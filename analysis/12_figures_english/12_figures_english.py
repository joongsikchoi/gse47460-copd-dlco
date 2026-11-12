"""Step 12: English-labelled manuscript figures for the JKDISS submission (tables/figures in English).
Runs the original Korean figure scripts (07_manuscript_figures.py, 11_figure6_deconvolution.py) after an ordered
string translation, so the plotted data and layout are identical; only labels change.
Usage: python 12_figures_english.py <07 script> <11 script> <analysis dir> <qc dir> <10_deconvolution dir> <out dir>"""
import sys, re, runpy
s07, s11, A, Q, D, out = sys.argv[1:7]
T07 = [("GSE47460 전체 표본 ", "All GSE47460 samples "), ("완전사례 주코호트 ", "Complete-case cohort "),
 ("PTX 표기 2명 포함 민감도 코호트 188 (126 · 62)", "Incl. 2 PTX entries: 188 (126 · 62)"),
 ("최종 주분석 코호트 177", "Final analysis cohort 177"), ("탐색 GPL14550 118 · 재현 GPL6480 59", "Discovery 118 · Replication 59"),
 ("다른 진단군 362명 제외", "Excluded: other diagnoses 362"), ("간질성폐질환 254 · 대조군 108", "ILD 254 · controls 108"),
 ("변수 불완전 34명 제외 (사유 중복)", "Excluded: incomplete data 34"), ("DLCO 결측 23 · 폐기종 결측 9", "DLCO 23 · emphysema 9 (overlap)"),
 ("흡연 결측 3 · FEV1 결측 2", "Smoking 3 · FEV1 2"), ("폐기종 문자 표기(PTX) 2", "Emphysema as text (PTX) 2"),
 ("배열 품질 이상치 9명 제외", "Excluded: array-quality outliers 9"), ("배열 간 상관 중앙값 < Q1 − 3×IQR", "Median array correlation < Q1 − 3×IQR"),
 ("탐색 GPL14550 (n = 118)", "Discovery GPL14550 (n = 118)"), ("재현 GPL6480 (n = 59)", "Replication GPL6480 (n = 59)"),
 (": FDR < 0.05 유전자 {sig.sum()}개", ": {sig.sum()} gene(s) at FDR < 0.05"), ("DLCO +10%p당 log2 발현 변화 (보정)", "Adjusted log2 expression change per +10 DLCO points"),
 ("전체 11,394개", "All 11,394 genes"), ("탐색 상위 100개", "Discovery top 100"), ("메타 FDR < 0.05 (7개)", "Meta FDR < 0.05 (7)"),
 ("플랫폼 간 t 통계량: 상위 100개 방향 일치 {conc*100:.0f}%", "Cross-platform t: top-100 direction concordance {conc*100:.0f}%"),
 ("탐색 t (GPL14550)", "Discovery t (GPL14550)"), ("재현 t (GPL6480)", "Replication t (GPL6480)"),
 ("탐색 내부 split-half", "Within-discovery split-half"), ("교차 플랫폼", "Cross-platform"), ("QC 전", "Before QC"), ("QC 후", "After QC"),
 ("상위 500개 유전자 효과 방향 일치율 (0.5 = 우연)", "Direction concordance of top 500 genes (0.5 = chance)"),
 ("재현 일치율: split-half는 중앙값과 5–95% 범위 (200회)", "Replication concordance: split-half median and 5–95% range (200 splits)"),
 ("AT2 (폐포)", "AT2"), ("폐포 모세혈관", "Alv. capillary"), ('"섬모"', '"Cilia"'), ("기도 분비", "Airway secretory"),
 ('"지방"', '"Adipose"'), ('"림프"', '"Lymphoid"'), ("조직 구성 표지 점수 (z 평균)", "Tissue marker score (mean z)"),
 ("QC 통과 177", "Passed QC 177"), ("QC 이상치 9", "QC outliers 9"),
 ("QC 이상치는 폐포 실질 표지가 낮고 비폐포 조직 표지가 높음", "QC outliers: low alveolar and high non-alveolar tissue markers"),
 ('"탐색")', '"discovery")'), ('"재현")', '"replication")'), ("{lab} 메타", "{lab} meta"),
 ("DLCO +10%p당 GSVA 점수 변화 (SD), FEV1·폐기종·나이·성별·흡연 보정, 95% CI", "GSVA score change (SD) per +10 DLCO points, adjusted, 95% CI"),
 ("H2: 사전 정의 폐포 세포유형 점수와 DLCO", "H2: prespecified alveolar cell-type scores and DLCO"),
 ("A: 임상변수", "A: clinical"), ("B: A + H2 점수 4개", "B: A + 4 H2 scores"), ("C: A + 탐색 100유전자 서명", "C: A + 100-gene signature"),
 ("재현 R² ", "Replication R² "), ("DLCO 설명력 R²", "DLCO R²"), ("탐색 학습 R²", "Discovery training R²"),
 ("재현 검증 R² (부트스트랩 95% CI)", "Replication R² (bootstrap 95% CI)"), (": 이상치 {(~ok).sum()}개", ": {(~ok).sum()} outliers"),
 ("다른 배열과의 상관 중앙값 (medcor)", "Median correlation with other arrays (medcor)"), ("PC1–5 제곱 Mahalanobis 거리", "Squared Mahalanobis distance, PC1–5"),
 ("'~')", "' to ')")]
T11 = [("a1.set_yticks(yt)", "a1.set_ylim(-11.8, 1.1); a1.set_yticks(yt)"), ("a2.set_yticks(yt)", "a2.set_ylim(-4.1, 0.5); a2.set_yticks(yt)"), ("폐포-모세혈관 합계", "Alveolar-capillary total"), ("aCap (aerocyte)", "aCap (aerocyte)"), ('"대식세포"', '"Macrophage"'),
 ('"기질세포"', '"Stromal"'), ('"림프구"', '"Lymphoid"'), ('"기도상피"', '"Airway epithelial"'), ('"메타 p"', '"Meta p"'),
 ("DLCO +10%p당 추정 세포 비율 변화 (%p), 95% CI", "Change in estimated fraction (points) per +10 DLCO points, 95% CI"),
 ("추정 세포 구성비와 잔여 DLCO", "Estimated cell fractions and residual DLCO"), ("탐색 (n = 118)", "Discovery (n = 118)"),
 ("재현 (n = 59)", "Replication (n = 59)"), ("역분산 가중 메타", "Inverse-variance meta"), ('"보정 전"', '"before"'),
 ('"구성비 보정 후"', '"after"'), ("보정 전 (주분석)", "Before adjustment (main)"), ("비폐포 구성비 보정 후", "After non-alveolar adjustment"),
 ("DLCO +10%p당 GSVA 점수 변화 (SD), 메타 95% CI", "GSVA score change (SD) per +10 DLCO points, meta 95% CI"),
 ("H2 점수의 구성비 보정 전후", "H2 scores before and after composition adjustment")]
def run(src, table, argv):
    code = open(src, encoding="utf8").read()
    for a, b in table:
        if a not in code: print("WARN missing:", a[:50])
        code = code.replace(a, b)
    left = [l for l in code.splitlines() if re.search("[가-힣]", l) and not l.lstrip().startswith("#") and '"""' not in l]
    if left: print("untranslated:", *left[:5], sep="\n  ")
    path = out + "/_tmp_" + src.split("/")[-1]; open(path, "w", encoding="utf8").write(code)
    sys.argv = [path] + argv; runpy.run_path(path, run_name="__main__")
run(s07, T07, [A, Q, out]); run(s11, T11, [D, out])
