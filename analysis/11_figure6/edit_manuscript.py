"""Integrate new analyses (deconvolution, robustness) into the existing manuscript docx (edit in place on a copy).
Usage: python edit_manuscript.py <unpacked dir> <gse dir> <figure6 png>"""
import sys, re, shutil, pandas as pd, numpy as np
U, G, FIG6 = sys.argv[1:4]
x = open(f'{U}/word/document.xml', encoding='utf8').read()

def esc(s): return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
def runs(text, rpr=''):
    out = []
    for i, part in enumerate(re.split(r'\[\[([^\]]+)\]\]', text)):
        if not part: continue
        if i % 2: out.append(f'<w:r><w:rPr><w:vertAlign w:val="superscript"/></w:rPr><w:t xml:space="preserve">{esc(part)}</w:t></w:r>')
        else: out.append(f'<w:r>{"<w:rPr>"+rpr+"</w:rPr>" if rpr else ""}<w:t xml:space="preserve">{esc(part)}</w:t></w:r>')
    return ''.join(out)
BODY = '<w:pPr><w:spacing w:after="120" w:line="340"/><w:ind w:firstLine="300"/><w:jc w:val="both"/></w:pPr>'
def body(text, lead=None):
    l = f'<w:r><w:rPr><w:b/><w:bCs/></w:rPr><w:t xml:space="preserve">{esc(lead)}</w:t></w:r>' if lead else ''
    return f'<w:p>{BODY}{l}{runs((" " if lead else "") + text)}</w:p>'
def abst(lead, text):
    return (f'<w:p><w:pPr><w:spacing w:after="100" w:line="320"/><w:jc w:val="left"/></w:pPr><w:r><w:rPr><w:b/><w:bCs/></w:rPr>'
            f'<w:t xml:space="preserve">{esc(lead)}</w:t></w:r>{runs(" " + text)}</w:p>')
def head(text, lvl): return f'<w:p><w:pPr><w:pStyle w:val="Heading{lvl}"/></w:pPr><w:r><w:t xml:space="preserve">{esc(text)}</w:t></w:r></w:p>'
def bullet(text): return (f'<w:p><w:pPr><w:pStyle w:val="ListParagraph"/><w:numPr><w:ilvl w:val="0"/><w:numId w:val="2"/></w:numPr>'
                          f'<w:spacing w:after="60" w:line="320"/></w:pPr>{runs(text)}</w:p>')
def ref(text): return (f'<w:p><w:pPr><w:pStyle w:val="ListParagraph"/><w:numPr><w:ilvl w:val="0"/><w:numId w:val="3"/></w:numPr>'
                       f'<w:spacing w:after="60" w:line="280"/></w:pPr><w:r><w:rPr><w:sz w:val="19"/><w:szCs w:val="19"/></w:rPr><w:t xml:space="preserve">{esc(text)}</w:t></w:r></w:p>')
def caption(label, text):
    return (f'<w:p><w:pPr><w:keepNext/><w:spacing w:after="80" w:before="200" w:line="280"/></w:pPr><w:r><w:rPr><w:b/><w:bCs/><w:sz w:val="19"/><w:szCs w:val="19"/></w:rPr>'
            f'<w:t xml:space="preserve">{esc(label)}</w:t></w:r><w:r><w:rPr><w:sz w:val="19"/><w:szCs w:val="19"/></w:rPr><w:t xml:space="preserve"> {esc(text)}</w:t></w:r></w:p>')
def note(text): return (f'<w:p><w:pPr><w:spacing w:after="220" w:before="50" w:line="250"/></w:pPr><w:r><w:rPr><w:color w:val="333333"/>'
                        f'<w:sz w:val="16"/><w:szCs w:val="16"/></w:rPr><w:t xml:space="preserve">{esc(text)}</w:t></w:r></w:p>')
NB = '<w:{0} w:val="none" w:color="FFFFFF" w:sz="0"/>'; SB = '<w:{0} w:val="single" w:color="000000" w:sz="4"/>'
def cell(t, w, top, bot, bold, jc):
    b = (SB if top else NB).format('top') + NB.format('left') + (SB if bot else NB).format('bottom') + NB.format('right')
    return (f'<w:tc><w:tcPr><w:tcW w:type="dxa" w:w="{w}"/><w:tcBorders>{b}</w:tcBorders><w:tcMar><w:top w:type="dxa" w:w="30"/><w:left w:type="dxa" w:w="70"/>'
            f'<w:bottom w:type="dxa" w:w="30"/><w:right w:type="dxa" w:w="70"/></w:tcMar></w:tcPr><w:p><w:pPr><w:spacing w:after="0" w:line="240"/><w:jc w:val="{jc}"/></w:pPr>'
            f'<w:r><w:rPr><w:rFonts w:ascii="Arial" w:cs="Arial" w:eastAsia="Malgun Gothic" w:hAnsi="Arial"/>{"<w:b/><w:bCs/>" if bold else ""}<w:sz w:val="15"/><w:szCs w:val="15"/></w:rPr>'
            f'<w:t xml:space="preserve">{esc(t)}</w:t></w:r></w:p></w:tc>')
def table(hdr, rows, widths, group_rows=()):
    assert sum(widths) == 9412, sum(widths)
    tb = ('<w:tbl><w:tblPr><w:tblW w:type="dxa" w:w="9412"/><w:tblBorders>' + ''.join(f'<w:{s} w:val="single" w:color="auto" w:sz="4"/>' for s in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV'])
          + '</w:tblBorders></w:tblPr><w:tblGrid>' + ''.join(f'<w:gridCol w:w="{w}"/>' for w in widths) + '</w:tblGrid>')
    jc = lambda j: 'left' if j == 0 else 'center'
    tb += '<w:tr><w:trPr><w:tblHeader/></w:trPr>' + ''.join(cell(h, w, 1, 1, 1, jc(j)) for j, (h, w) in enumerate(zip(hdr, widths))) + '</w:tr>'
    for i, r in enumerate(rows):
        last = i == len(rows) - 1; g = i in group_rows
        tb += '<w:tr><w:trPr><w:cantSplit/></w:trPr>' + ''.join(cell(str(c), w, 0, last, g, jc(j)) for j, (c, w) in enumerate(zip(r, widths))) + '</w:tr>'
    return tb + '</w:tbl><w:p><w:pPr><w:spacing w:after="120"/></w:pPr></w:p>'

def pbounds(key, start=0):
    i = x.find(key, start); assert i >= 0 and x.count(key) == 1, (key, x.count(key))
    s = max(x.rfind('<w:p>', 0, i), x.rfind('<w:p ', 0, i)); e = x.find('</w:p>', i) + 6; return s, e
def rep_para(key, new):
    global x; s, e = pbounds(key); x = x[:s] + new + x[e:]
def after_para(key, new):
    global x; s, e = pbounds(key); x = x[:e] + new + x[e:]
def before_para(key, new):
    global x; s, e = pbounds(key); x = x[:s] + new + x[s:]
def rep(old, new):
    global x; assert x.count(old) == 1, (old, x.count(old)); x = x.replace(old, new)

# ---------- data ----------
pr = pd.read_csv(f'{G}/out10/deconv_proportions_vs_DLCO.csv').set_index('proportion')
dc = pd.read_csv(f'{G}/out10/deconv_proportions_186.csv')
val = pd.read_csv(f'{G}/out10/deconv_validation_vs_gsva.csv')
h2d = pd.read_csv(f'{G}/out10/H5_deconv_H2.csv').set_index('set')
mv = pd.read_csv(f'{G}/out9/A_model_variants.csv').set_index('model')
mz = pd.read_csv(f'{G}/out9/B_meanz_meta.csv').set_index('set')
loo = pd.read_csv(f'{G}/out9/C_loo_H2_meta.csv').set_index('set')
cl = pd.read_csv(f'{G}/out9/D_lowtertile_classifier.csv')
cam = pd.read_csv(f'{G}/out9/F_camera_celltype_sets.csv'); cam0 = pd.read_csv(f'{G}/out9/F_camera_celltype_sets_fixedcor0.01.csv')
sens = pd.read_csv(f'{G}/out4b/H2_meta_sensitivity.csv')
f2 = lambda v: f'{v:.2f}'; f3 = lambda v: f'{v:.3f}'
def fp(p): return '< 0.001' if p < 0.001 else (f'{p:.4f}' if p < 0.01 else f'{p:.3f}')

# ---------- cover ----------
rep('본문, 표 4개, 그림 5개, 부록(방법 S1, 표 S1–S9, 그림 S1–S5)', '본문, 표 5개, 그림 6개, 부록(방법 S1–S2, 표 S1–S11, 그림 S1–S5)')
rep('2026년 10월 21일', '2026년 10월 29일 (통합 수정본)')

# ---------- abstract ----------
rep('사전에 정한 배열 품질 기준으로 9명을 제외한', '1차 분석 후 도입한 배열 품질 기준(배열 간 상관 극단 이상치)으로 9명을 제외한')
rep('세포 구성 대리변수 보정 후 유지 여부(H5)를 검정했다.', '세포 구성 보정 후 유지 여부(H5)를 검정했다. 세포 구성은 같은 연구진의 대조군·COPD 폐 단일세포 자료(GSE136831)를 참조로 한 표지 기반 디컨볼루션으로 추정했다.')
rep('품질 기준 도입 전 탐색 코호트에서 관찰된 375개 유의 유전자는 폐포 실질 비율이 낮은 배열 7개가 만든 신호였다.',
    '디컨볼루션으로 추정한 폐포-모세혈관 세포(AT1, AT2, aCap, gCap) 비율은 DLCO +10%p당 0.96%p 높았고(메타 p = 0.0013), aerocyte(aCap) 비율의 연관이 가장 일관되었다(메타 p = 0.0006). 비폐포 세포 비율을 보정하자 H2 점수의 연관은 약해졌다(AT1 메타 p 0.0026 → 0.076). 품질 기준 도입 전 탐색 코호트의 375개 유의 유전자는 폐포 세포 비율이 낮은 배열 7개가 만든 신호였다.')
rep_para('결론:', abst('결론:', 'COPD 수술 폐조직에서 기류제한과 폐기종으로 설명되지 않는 DLCO 변동은 특정 유전자의 발현 조절로는 재현성 있게 설명되지 않았다. 대신 같은 FEV1과 폐기종에서 DLCO가 낮은 환자의 조직에는 가스교환을 담당하는 폐포 모세혈관(특히 aerocyte)과 1형 폐포상피 세포가 적었다. 잔여 확산능 저하는 CT로 보이지 않는 폐포-모세혈관 단위의 소실을 반영할 수 있으며, 벌크 조직 연구에서 세포 구성을 고려하지 않으면 거짓 양성 유전자 신호가 생길 수 있다.'))
rep('일산화탄소 확산능, 유전자 발현, 폐포 모세혈관, 재현성', '일산화탄소 확산능, 유전자 발현, 폐포 모세혈관, 세포 구성 디컨볼루션, 재현성')

# ---------- methods ----------
after_para('벌크 조직의 발현 차이가',
    head('단일세포 참조 디컨볼루션 (추가 분석)', 3) +
    body('유전자군 점수는 세포 비율의 간접 지표이므로, 사전 분석 계획에 없던 추가 분석으로 단일세포 참조 기반 세포 비율을 추정했다. 참조는 같은 연구진의 폐 단일세포 RNA-seq 자료 GSE136831에서 대조군과 COPD 세포 163,246개(46명)를 원 저자의 세포유형 주석에 따라 17개 세포군으로 묶고, 세포별 계수를 총 UMI로 나눈 값의 세포군 평균(CP10K)으로 만들었다[[19]]. 폐 모세혈관은 원 자료의 Capillary A(aerocyte, aCap)와 Capillary B(gCap) 구분을 따랐다. 세포군마다 다른 세포군 최댓값 대비 log2 특이도가 1을 넘는 표지 유전자를 배열에서 최대 20개(세포군당 7–20개) 골라, Hunt 등의 표지 기반 추정식(dtangle, γ = 1)으로 표본별 상대 비율을 구했다[[20]](부록 방법 S2). 추정 비율(%)을 H2 점수와 같은 공변량 모형의 종속변수로 두어 DLCO 계수를 플랫폼별로 추정하고 역분산 가중으로 결합했다. 또한 디컨볼루션으로 추정한 비폐포 세포 비율(대식세포, 림프구, 기도상피, 기질세포)을 공변량으로 추가해 H2 점수와 유전자 수준 결과를 다시 적합했다(H5 재검정).') +
    head('추가 강건성 분석', 3) +
    body('H2 결과의 안정성을 보기 위해 (1) 표본 1개씩 제외한 메타분석(177회), (2) GSVA 대신 표지 유전자 z 점수 평균, (3) 나이·성별 보정 제외 및 채취 엽(상엽 여부) 추가 모형, (4) 유전자 간 상관을 추정한 camera 경쟁 유전자군 검정[[21]]과 GO 생물학적 과정 3,811개에 대한 camera 검정, (5) 잔여 DLCO 하위 1/3 여부의 분류 성능(탐색 학습, 재현 AUC)을 평가했다.'))
rep('통계 분석은 R(limma, GSVA, sva)과 Python(pandas, numpy, scipy, matplotlib)으로 수행했다.',
    '통계 분석은 R(limma, GSVA, sva)과 Python(pandas, numpy, scipy, matplotlib)으로 수행했고, 단일세포 참조 집계와 디컨볼루션은 JavaScript와 Python으로 수행했다.')

# ---------- results ----------
qc = dc.assign(Alv=dc[['AT1', 'AT2', 'CapA', 'CapB']].sum(axis=1))
m14 = qc[(qc.platform == 'GPL14550')].groupby('qc_outlier').Alv.median() * 100
after_para('Agilent 자체 스캔 품질 지표는', body(
    f'단일세포 참조 디컨볼루션도 같은 결론을 보였다. GPL14550의 이상치 7개 배열은 추정 폐포-모세혈관 세포 비율 중앙값이 {m14[True]:.1f}%로 나머지 118개({m14[False]:.1f}%)의 약 1/3이었다(Mann–Whitney p < 0.001; 부록 표 S10). GPL6480의 이상치 2개 중 GSM1149972는 5.7%로 낮았으나 GSM1150016은 42.5%로 오히려 높아, 이 배열은 조직 구성 외의 이유로 다른 배열과의 상관이 낮았던 것으로 보인다.'))

r = h2d
after_para('비폐포 세포 4종 점수를 추가 보정해도', body(
    '단일세포 참조 디컨볼루션 추정치는 GSVA 점수와 잘 맞았다(Spearman ρ: AT1 0.71–0.78, AT2 0.75–0.85, 모세혈관 대 gCap 0.84–0.87, 대식세포 0.78–0.85; 부록 표 S10 주석). '
    f'디컨볼루션으로 추정한 비폐포 세포 비율을 보정하자 H2 점수의 연관은 크게 약해졌다. 메타 β는 AT1 {r.loc["AT1","meta_beta"]:.3f} → {r.loc["AT1","meta_beta_adj"]:.3f} SD(p {fp(r.loc["AT1","meta_p"])} → {fp(r.loc["AT1","meta_p_adj"])}), '
    f'gCap {r.loc["gCap","meta_beta"]:.3f} → {r.loc["gCap","meta_beta_adj"]:.3f} SD(p {fp(r.loc["gCap","meta_p"])} → {fp(r.loc["gCap","meta_p_adj"])}), aCap {r.loc["aCap","meta_beta"]:.3f} → {r.loc["aCap","meta_beta_adj"]:.3f} SD(p {fp(r.loc["aCap","meta_p"])} → {fp(r.loc["aCap","meta_p_adj"])})였고, '
    '효과 크기는 탐색 코호트에서 원래의 23–40%, 재현 코호트에서 48–74%만 남았다(그림 6B). 유전자 수준에서는 탐색 상위 100개 유전자의 방향이 탐색 100%, 재현 83% 유지되었고 효과 크기 중앙값은 원래의 76%와 69%였으며, 보정 후 FDR < 0.05 유전자는 탐색 1개, 재현 0개였다. 따라서 H5는 유전자 수준에서는 대체로 유지되었으나, H2 점수의 연관은 대부분 세포 구성 차이로 설명되었다.') +
    head('3.7 추정 세포 구성비와 잔여 DLCO (추가 분석)', 2) +
    body('H2 점수의 감쇠가 세포 양의 차이를 뜻하는지 직접 확인하기 위해 추정 세포 비율 자체를 잔여 DLCO와 연관지었다(표 5, 그림 6A). FEV1, 폐기종, 나이, 성별, 흡연을 보정한 뒤 DLCO가 10%p 높을수록 폐포-모세혈관 세포(AT1 + AT2 + aCap + gCap) 비율은 탐색 1.02%p(p = 0.003), 재현 0.73%p(p = 0.27) 높았고, 메타 추정치는 0.96%p(p = 0.0013, 이질성 p = 0.70)였다. 개별 세포군 중에서는 aerocyte(aCap) 비율(메타 0.17%p, p = 0.0006; 탐색 p = 0.004, 재현 p = 0.09)과 AT1 비율(메타 0.30%p, p = 0.0025)의 연관이 가장 일관되었고, 두 코호트의 추정치가 거의 같았다. 반대로 대식세포 비율은 DLCO가 높을수록 낮았다(메타 −0.68%p, p = 0.033). AT2, 림프구, 기도상피, 기질세포 비율은 연관이 없었다. 즉 같은 FEV1과 폐기종에서 DLCO가 낮은 환자의 조직은 가스교환 단위를 이루는 세포가 적고 대식세포가 상대적으로 많은 구성을 보였다.') +
    caption('표 5.', '단일세포 참조 디컨볼루션으로 추정한 세포 비율과 DLCO의 보정 연관 (추가 분석, n = 177)') +
    table(['세포군', '비율 중앙값 % (탐색/재현)', '탐색 β (p)', '재현 β (p)', '메타 β', '메타 p', '이질성 p'],
          [[lab, f'{pr.loc[k,"disc_median_pct"]:.1f} / {pr.loc[k,"rep_median_pct"]:.1f}', f'{pr.loc[k,"disc_beta_pct_per10"]:.2f} ({fp(pr.loc[k,"disc_p"])})',
            f'{pr.loc[k,"rep_beta_pct_per10"]:.2f} ({fp(pr.loc[k,"rep_p"])})', f'{pr.loc[k,"meta_beta"]:.2f}', fp(pr.loc[k,"meta_p"]), f'{pr.loc[k,"Q_p"]:.2f}']
           for k, lab in [('AlvCap', '폐포-모세혈관 합계'), ('AT1', 'AT1'), ('AT2', 'AT2'), ('CapA', 'aCap (aerocyte)'), ('CapB', 'gCap'), ('Cap', '모세혈관 합계'),
                          ('Macrophage', '대식세포'), ('Stromal', '기질세포'), ('Lymphoid', '림프구'), ('Airway', '기도상피')]],
          [1880, 1560, 1420, 1420, 1000, 1000, 1132]) +
    note('β: DLCO 10%p 증가당 추정 세포 비율 변화(%p), FEV1·√폐기종·나이·성별·현재 흡연 보정. 메타: 역분산 가중 고정효과, 이질성: Cochran Q. 비율은 GSE136831 대조군·COPD 세포 참조와 dtangle 방식으로 추정한 상대 비율로, 절대값보다 표본 간 차이를 해석한다. 기질세포 = 섬유아세포 + 평활근·혈관주위세포, 림프구 = T/NK + B/형질세포, 기도상피 = 기저·클럽·배상·섬모세포 등. 사전 계획에 없던 추가 분석이다.') +
    f'<w:p><w:pPr><w:keepNext/><w:spacing w:after="60" w:before="160"/><w:jc w:val="center"/></w:pPr><w:r><w:drawing><wp:inline distT="0" distB="0" distL="0" distR="0"><wp:extent cx="5900000" cy="2675000"/><wp:effectExtent t="0" r="0" b="0" l="0"/><wp:docPr id="40" name="Figure6" descr="" title=""/><wp:cNvGraphicFramePr><a:graphicFrameLocks xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" noChangeAspect="1"/></wp:cNvGraphicFramePr><a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:nvPicPr><pic:cNvPr id="0" name="Figure6.png" descr=""/><pic:cNvPicPr><a:picLocks noChangeAspect="1" noChangeArrowheads="1"/></pic:cNvPicPr></pic:nvPicPr><pic:blipFill><a:blip r:embed="rId40" cstate="none"/><a:srcRect/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr bwMode="auto"><a:xfrm><a:off x="0" y="0"/><a:ext cx="5900000" cy="2675000"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>' +
    caption('그림 6.', '단일세포 참조 디컨볼루션 결과 (추가 분석). (A) 추정 세포 비율과 DLCO의 보정 연관(FEV1·폐기종·나이·성별·흡연 보정, 95% CI). 파랑: 탐색, 주황: 재현, 검은 마름모: 역분산 가중 메타분석. (B) H2 점수의 메타 추정치: 주분석(회색)과 비폐포 세포 비율(대식세포, 림프구, 기도상피, 기질세포) 보정 후(검정).'))

rep('3.7 예측 설명력 (H4)', '3.8 예측 설명력 (H4)')
rep('3.8 배치와 민감도 분석', '3.9 배치와 민감도 분석')
c = cl.set_index('model').loc['4 H2 scores']
after_para('재현 코호트에서 임상변수 모형 A', body(f'잔여 DLCO 하위 1/3 여부를 분류하는 보조 분석에서도 H2 점수 4개의 재현 코호트 AUC는 {c.AUC_replication:.2f}(95% CI {c.CI_low:.2f}–{c.CI_high:.2f})로 우연 수준이었다(부록 표 S11).'))
after_para('스캔 날짜(GPL14550', body(
    f'H2 메타분석은 표본 1개 제외 분석에서 AT1과 gCap은 177회 모두 p < 0.05였고(최대 p {loo.loc["AT1","meta_p_max"]:.3f}, {loo.loc["gCap","meta_p_max"]:.3f}), aCap은 {loo.loc["aCap","frac_p_below_05"]*100:.1f}%에서 유지되었다. GSVA 대신 표지 z 점수 평균을 써도 AT1(p = {mz.loc["AT1","meta_p"]:.3f}), aCap(p = {mz.loc["aCap","meta_p"]:.3f}), gCap(p = {mz.loc["gCap","meta_p"]:.3f})의 양의 연관이 유지되었고, 나이·성별 보정 제외나 상엽 여부 추가 모형에서도 같았다(부록 표 S11). '
    '반면 유전자 간 상관을 추정한 camera 경쟁 검정에서는 어떤 H2 유전자군도 유의하지 않았고(탐색 p 0.06–0.72), GO 생물학적 과정 3,811개 중 두 코호트에서 모두 FDR < 0.05인 경로는 없었다. 즉 폐포 세포 신호는 표본 간 세포 양의 차이로서는 안정적이었으나, 해당 유전자군이 다른 유전자보다 특별히 조절되는 신호는 아니었다.'))

# Table 4: update H5 row, add a row for deconvolution
rep('세포 구성 점수 추가 보정', '세포 구성 보정 (유전자군 점수, 단일세포 디컨볼루션)')
rep('상위 100개 방향 100% / 87% 유지; H2 점수는 탐색에서 약 50% 감쇠', '유전자: 상위 100개 방향 100% / 83% 유지; H2 점수: 디컨볼루션 보정 시 AT1 메타 p 0.0026 → 0.076')
rep('부분 유지', '유전자 유지, H2 점수 감쇠')
i = x.find('가설 검정 요약'); t = x.find('<w:tbl>', i); te = x.find('</w:tbl>', t)
rws = [m for m in re.finditer(r'<w:tr>.*?</w:tr>', x[t:te], re.S)]
last = rws[-1].group(0)
newrow = last
for o, n in [('>H5<', '>추가<'), ('세포 구성 보정 (유전자군 점수, 단일세포 디컨볼루션)', '추정 세포 비율 회귀'),
             ('유전자: 상위 100개 방향 100% / 83% 유지; H2 점수: 디컨볼루션 보정 시 AT1 메타 p 0.0026 → 0.076', '폐포-모세혈관 비율 메타 +0.96%p/10%p(p = 0.0013); aCap p = 0.0006, AT1 p = 0.0025'),
             ('유전자 유지, H2 점수 감쇠', '지지 (탐색적)')]:
    assert o in newrow, o; newrow = newrow.replace(o, n)
crit = re.findall(r'<w:t xml:space="preserve">([^<]*)</w:t>', last)[2]
newrow = newrow.replace(f'>{crit}<', '>p &lt; 0.05 (사후 분석)<', 1)
oldlast = last.replace('<w:bottom w:val="single" w:color="000000" w:sz="4"/>', '<w:bottom w:val="none" w:color="FFFFFF" w:sz="0"/>')
a, b = t + rws[-1].start(), t + rws[-1].end(); x = x[:a] + oldlast + newrow + x[b:]

# ---------- discussion ----------
rep_para('본 연구는 COPD 수술 폐조직에서', body('본 연구는 COPD 수술 폐조직에서 기류제한과 CT 폐기종으로 설명되지 않는 DLCO 변동과 연관된 발현 특성을 두 독립 플랫폼에서 검토했다. 개별 유전자 수준에서는 재현 가능한 신호가 없었으나, 사전 정의한 1형 폐포상피와 폐포 모세혈관 세포 점수는 두 코호트에서 같은 방향과 비슷한 크기로 잔여 DLCO와 양의 연관을 보였다. 단일세포 참조로 세포 비율을 직접 추정하자 이 연관은 세포 양의 차이로 정리되었다. 같은 FEV1과 폐기종에서 DLCO가 낮은 환자의 조직에는 폐포-모세혈관 세포, 특히 aerocyte와 AT1 세포가 적었고, 이 구성비를 보정하면 유전자군 점수의 연관은 대부분 사라졌다. 즉 잔여 DLCO에 대응하는 조직 특성은 특정 유전자의 발현 조절이 아니라 가스교환 단위를 이루는 세포의 양이며, 메타분석에서 가장 강한 개별 유전자가 gCap 표지자 GPIHBP1이었던 것도 같은 맥락이다.'))
rep('폐포 모세혈관·AT1 표지가 적다는', '폐포 모세혈관(특히 가스교환을 맡는 aerocyte)과 AT1 세포가 적다는')
rep_para('다만 이 신호가', body('벌크 자료만으로는 이 신호가 세포 수의 차이인지 세포당 발현의 변화인지 구분하기 어렵지만, 디컨볼루션 결과는 전자를 지지한다. 추정 세포 비율 자체가 잔여 DLCO와 연관되었고, 비폐포 세포 비율을 보정하자 H2 점수의 연관이 탐색 코호트에서 원래의 23–40%로 줄었다. 유전자 간 상관을 고려한 camera 경쟁 검정에서 H2 유전자군이 다른 유전자보다 특별히 강하게 연관되지 않은 것도, 신호가 특정 경로의 조절보다 조직 내 세포 양의 차이에서 온다는 해석과 맞는다. 폐포 실질이 적은 조직은 실제 폐의 미세 구조 손상을 반영할 수도 있고 채취 위치의 차이일 수도 있으나, 추정 비율이 FEV1과 폐기종을 보정한 DLCO와 연관되었다는 점은 적어도 일부가 환자의 폐 상태를 반영함을 시사한다. 대식세포 비율이 DLCO와 음의 연관을 보인 것은 폐포 내 염증세포 축적을 반영할 수 있다. AT2 점수와 비율이 연관을 보이지 않은 것은 COPD 폐에서 AT2 아형의 구성 변화가 크다는 단일세포 결과[[4]]와 함께, 전체 AT2 합산으로는 아형별 반대 방향 변화가 상쇄될 수 있음을 시사한다.'))
rep('이 배열들은 제조사 스캔 품질 지표로는 걸러지지 않았다', '탐색 코호트의 이 배열들은 추정 폐포-모세혈관 세포 비율이 나머지의 약 1/3이었지만 제조사 스캔 품질 지표로는 걸러지지 않았다')
rep_para('세포 구성은 참조 기반', bullet('세포 구성은 단일세포 참조(GSE136831) 기반 표지 디컨볼루션으로 추정했으나, 단일세포와 마이크로어레이의 측정 차이로 절대 비율은 근사치다. 디컨볼루션 분석은 사전 계획에 없던 추가 분석이며, 조직 형태계측 등 독립적 측정으로 확인할 필요가 있다.'))
rep_para('COPD 수술 폐조직에서 FEV1과 CT 폐기종으로', body('COPD 수술 폐조직에서 FEV1과 CT 폐기종으로 설명되지 않는 DLCO 변동은 특정 유전자의 발현 조절로는 재현성 있게 설명되지 않았다. 대신 DLCO가 낮은 환자의 조직에는 가스교환을 담당하는 폐포 모세혈관(특히 aerocyte)과 1형 폐포상피 세포가 적었고, 이 관계는 두 플랫폼에서 같은 방향이었다. 잔여 확산능 저하는 CT로 보이지 않는 폐포-모세혈관 단위의 소실을 반영할 가능성이 있다. 또한 벌크 폐조직 발현 연구에서는 조직 구성 이상치가 수백 개의 거짓 연관을 만들 수 있으므로, 독립 코호트 재현과 세포 구성 평가가 필수적이다.'))
s, e = pbounds('R 6개, Python')
x = x[:s] + body('분석에 사용한 자료는 NCBI GEO(GSE47460, 단일세포 참조 GSE136831)에서 공개되어 있다. 코호트 구성, limma·GSVA 분석, 재현 진단, 배열 품질 지표, 표·그림 생성 코드와 실행 순서는 연구 자료 폴더의 code/ 디렉터리(README.md)에, 디컨볼루션과 추가 강건성 분석 코드와 결과는 analysis/08_figures–11_figure6 디렉터리에 있다.') + x[e:]

# ---------- references ----------
rep('[저자 확인 필요].', 'Zheng W, Zou M, Hu X, et al.')
rep('Commun Biol. 2025.', 'Commun Biol. 2025;8:93.')
after_para('Whitlock MC', ref('Adams TS, Schupp JC, Poli S, et al. Single-cell RNA-seq reveals ectopic and aberrant lung-resident cell populations in idiopathic pulmonary fibrosis. Sci Adv. 2020;6(28):eaba1983. doi:10.1126/sciadv.aba1983') +
           ref('Hunt GJ, Freytag S, Bahlo M, Gagnon-Bartsch JA. dtangle: accurate and robust cell type deconvolution. Bioinformatics. 2019;35(12):2093-2099. doi:10.1093/bioinformatics/bty926') +
           ref('Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Res. 2012;40(17):e133.'))

# ---------- supplement ----------
rep('부록 방법 S1 · 부록 표 S1–S9 · 부록 그림 S1–S5', '부록 방법 S1–S2 · 부록 표 S1–S11 · 부록 그림 S1–S5')
before_para('부록 표 S1.', head('부록 방법 S2. 단일세포 참조 디컨볼루션', 2) +
    body('GSE136831의 원시 계수 행렬(45,947개 유전자 × 312,928개 세포)과 세포 메타데이터를 읽어, 특발성폐섬유증 세포를 제외한 대조군·COPD 세포를 원 저자의 세포유형 주석으로 18개 군으로 묶고 폐포·간질 대식세포를 합쳐 17개 군을 만들었다. 세포마다 계수를 총 UMI로 나누고 세포군별로 평균하여 10⁴를 곱했다(CP10K). 각 세포군에서 CP10K ≥ 0.5이고 다른 군 최댓값 대비 log2(x + 1) 차이가 0.5를 넘는 유전자 상위 100개를 모아 참조 행렬을 만들었다. GEO는 분석 서버에서 접근할 수 없어, 2.0 GB 계수 파일을 연구자 PC의 브라우저에서 분할 요청(Range)으로 스트리밍 집계했다(analysis/10_deconvolution/README_reference_aggregation.md).', lead='참조 구축.') +
    body('플랫폼별로 같은 심볼의 프로브를 평균하고, 참조와 배열에 공통인 유전자 중 log2 특이도 > 1, 참조 CP10K ≥ 1인 표지를 세포군당 최대 20개 골랐다. 표본 j, 세포군 k의 지수 a_jk는 표지 유전자의 (y_gj − log2(r_gk + 1)) 평균이며, 비율은 p_jk = 2^(a_jk/γ) / Σ_l 2^(a_jl/γ), γ = 1로 계산했다(민감도 분석 γ = 0.7). 프로브별 신호 차이 때문에 절대 비율은 근사치이며, 모형에는 플랫폼 내 표본 간 변이가 들어간다. 비음수 최소제곱 디컨볼루션은 모세혈관 비율이 0으로 수렴하는 등 적합이 불량하여 사용하지 않았다.', lead='비율 추정.'))

# S10
dd = dc.copy(); dd['Alv'] = dd[['AT1', 'AT2', 'CapA', 'CapB']].sum(axis=1); dd['Lym'] = dd.T_NK + dd.B_Plasma; dd['Str'] = dd.Fibroblast + dd.SMC_Pericyte
rows10 = []
for _, q in dd[dd.qc_outlier].iterrows():
    rows10.append([q.sample_id, q.platform.replace('GPL', ''), f'{q.Alv*100:.1f}', f'{q.AT1*100:.1f}', f'{(q.CapA+q.CapB)*100:.1f}', f'{q.Macrophage*100:.1f}', f'{q.Str*100:.1f}', f'{q.Lym*100:.1f}'])
for p in ['GPL14550', 'GPL6480']:
    s = dd[(dd.platform == p) & ~dd.qc_outlier]
    rows10.append([f'유지 배열 중앙값 (n = {len(s)})', p.replace('GPL', ''), *[f'{s[c].median()*100:.1f}' for c in ['Alv', 'AT1']], f'{(s.CapA+s.CapB).median()*100:.1f}', *[f'{s[c].median()*100:.1f}' for c in ['Macrophage', 'Str', 'Lym']]])
vtxt = '; '.join(f'{r.platform}: AT1 {r.AT1_vs_GSVA_AT1:.2f}, AT2 {r.AT2_vs_GSVA_AT2:.2f}, 모세혈관–gCap {r.Cap_vs_GSVA_gCap:.2f}, 모세혈관–aCap {r.Cap_vs_GSVA_aCap:.2f}, 대식세포 {r.Mac_vs_GSVA_Mac:.2f}' for r in val.itertuples())
# S11
def cellb(b, p): return f'{b:.3f} ({fp(p)})'
S = sens.set_index(['cohort', 'set'])
sets = ['AT1', 'gCap', 'aCap', 'AT2']
rows11 = [['주분석 (GSVA, 메타)'] + [cellb(S.loc[('primary_QC', k), 'meta_beta_per10'], S.loc[('primary_QC', k), 'meta_p']) for k in sets],
          ['QC + GOLD 0 제외 (115/55)'] + [cellb(S.loc[('QC_GOLD1to4', k), 'meta_beta_per10'], S.loc[('QC_GOLD1to4', k), 'meta_p']) for k in sets],
          ['QC + 성별 불일치 제외 (115/58)'] + [cellb(S.loc[('QC_sexconcordant', k), 'meta_beta_per10'], S.loc[('QC_sexconcordant', k), 'meta_p']) for k in sets],
          ['품질관리 전 186명'] + [cellb(S.loc[('noQC_186', k), 'meta_beta_per10'], S.loc[('noQC_186', k), 'meta_p']) for k in sets],
          ['나이·성별 보정 제외'] + [cellb(mv.loc['no_age_sex', f'{k}_meta_beta'], mv.loc['no_age_sex', f'{k}_meta_p']) for k in sets],
          ['+ 상엽 여부 보정'] + [cellb(mv.loc['plus_upper_lobe', f'{k}_meta_beta'], mv.loc['plus_upper_lobe', f'{k}_meta_p']) for k in sets],
          ['표지 z 점수 평균'] + [cellb(mz.loc[k, 'meta'], mz.loc[k, 'meta_p']) for k in sets],
          ['표본 1개 제외: 최대 p (p < 0.05 비율)'] + [f'{loo.loc[k,"meta_p_max"]:.3f} ({loo.loc[k,"frac_p_below_05"]*100:.1f}%)' for k in ['AT1', 'gCap', 'aCap']] + ['—'],
          ['디컨볼루션 비폐포 비율 보정'] + [cellb(h2d.loc[k, 'meta_beta_adj'], h2d.loc[k, 'meta_p_adj']) for k in sets],
          ['camera, 상관 추정: 탐색 / 재현 p'] + [f'{cam.set_index(["platform","set"]).loc[("GPL14550","H2_"+k),"PValue"]:.3f} / {cam.set_index(["platform","set"]).loc[("GPL6480","H2_"+k),"PValue"]:.3f}' for k in sets],
          ['camera, 상관 0.01 고정: 탐색 / 재현 p'] + [f'{fp(cam0.set_index(["platform","set"]).loc[("GPL14550","H2_"+k),"PValue"])} / {fp(cam0.set_index(["platform","set"]).loc[("GPL6480","H2_"+k),"PValue"])}' for k in sets]]
cl2 = cl.set_index('model')
after_para('H5 구성 보정용 비폐포',
    caption('부록 표 S10.', '품질관리 제외 배열의 디컨볼루션 추정 세포 비율 (%)') +
    table(['GSM', 'GPL', '폐포-모세혈관', 'AT1', '모세혈관', '대식세포', '기질세포', '림프구'], rows10, [2112, 800, 1300, 900, 1100, 1100, 1100, 1000]) +
    note(f'폐포-모세혈관 = AT1 + AT2 + aCap + gCap, 모세혈관 = aCap + gCap. GPL14550 이상치 7개 대 유지 118개의 폐포-모세혈관 비율 Mann–Whitney p < 0.001. 디컨볼루션 추정치와 GSVA 점수의 Spearman 상관 — {vtxt}.') +
    caption('부록 표 S11.', 'H2 유전자군 결과의 강건성 (QC 177명, DLCO +10%p당 메타 β (p))') +
    table(['분석', 'AT1', 'gCap', 'aCap', 'AT2'], rows11, [2812, 1650, 1650, 1650, 1650]) +
    note(f'β: GSVA(또는 z 평균) 점수의 SD 변화. camera는 경쟁 유전자군 검정으로, 상관 0.01 고정(기본값)은 유전자 간 상관을 과소평가해 p가 작아질 수 있어 상관 추정 결과를 주결과로 본다. GO 생물학적 과정 3,811개 camera(상관 추정)에서 탐색·재현 모두 FDR < 0.05인 경로는 0개였다. 잔여 DLCO 하위 1/3 분류(탐색 학습): H2 점수 4개 재현 AUC {cl2.loc["4 H2 scores","AUC_replication"]:.2f}(95% CI {cl2.loc["4 H2 scores","CI_low"]:.2f}–{cl2.loc["4 H2 scores","CI_high"]:.2f}), AT1 점수 단독 {cl2.loc["AT1 score only","AUC_replication"]:.2f}({cl2.loc["AT1 score only","CI_low"]:.2f}–{cl2.loc["AT1 score only","CI_high"]:.2f}).'))

open(f'{U}/word/document.xml', 'w', encoding='utf8').write(x)
rel = open(f'{U}/word/_rels/document.xml.rels', encoding='utf8').read()
if 'rId40' not in rel:
    rel = rel.replace('</Relationships>', '<Relationship Id="rId40" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/figure6_deconvolution.png"/></Relationships>')
    open(f'{U}/word/_rels/document.xml.rels', 'w', encoding='utf8').write(rel)
shutil.copy(FIG6, f'{U}/word/media/figure6_deconvolution.png')
print('done')
