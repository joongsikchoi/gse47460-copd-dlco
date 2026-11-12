"""
GSE47460 COPD DLCO study - Step 7: manuscript figures (Korean labels) from saved analysis outputs.

Main figures
  Figure1_cohort_flow.png        cohort selection 582 -> 177
  Figure2_H1_volcano_concordance cohort QC 177: discovery / replication volcano + cross-platform t plot
  Figure3_replication_QC.png     cross-platform vs split-half concordance before/after QC + tissue scores
  Figure4_H2_forest.png          H2 cell-type scores: discovery, replication, inverse-variance meta
  Figure5_H4_prediction.png      out-of-sample R2 in the replication platform
Supplementary
  FigureS5_array_QC_metrics.png  medcor vs PC_maha with the 9 excluded arrays

Usage: python 07_manuscript_figures.py <analysis dir> <qc_array dir> <out dir>
  <analysis dir> = folder holding 01_cohort/, 03_diagnostics/, 04_main/, 05_raw_qc/ outputs
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

A, Q, out = map(Path, sys.argv[1:4])
OUTL = {"GSM1150147","GSM1150187","GSM1150199","GSM1150220","GSM1150226","GSM1150251","GSM1150465","GSM1149972","GSM1150016"}

out.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, INK, MUTED, GRID, GREY = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6b68", "#e6e6e3", "#c9c9c4"
plt.rcParams.update({
    "font.family": ["Malgun Gothic", "Noto Sans CJK KR", "AppleGothic", "NanumGothic", "DejaVu Sans"],
    "axes.unicode_minus": False, "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight"})


def tag(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


# ---------------- Figure 1: cohort flow ----------------
fl = pd.read_csv(A / "01_cohort/cohort_flow.csv").set_index("Step")
fig, ax = plt.subplots(figsize=(7.2, 6.4))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")


def box(x, y, w, h, lines, main=False, side=False):
    fc = "#e8f0fb" if main else ("#f3f3f1" if side else "white")
    ec = BLUE if main else MUTED
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=1.6 if main else 0.9))
    for i, (t, bold) in enumerate(lines):
        ax.text(x + (1.5 if side else w / 2), y + h - 3.2 - i * 3.6, t, ha="left" if side else "center", va="center",
                fontsize=9 if bold else 8, fontweight="bold" if bold else "normal", color=INK if bold else MUTED)


r = lambda k: fl.loc[k]
box(5, 88, 42, 9, [(f"GSE47460 전체 표본 {r('All GSE47460 samples').Total}", True),
                   (f"GPL14550 {r('All GSE47460 samples').GPL14550} · GPL6480 {r('All GSE47460 samples').GPL6480}", False)])
box(5, 66, 42, 9, [(f"COPD {r('COPD').Total}", True), (f"GPL14550 {r('COPD').GPL14550} · GPL6480 {r('COPD').GPL6480}", False)])
box(5, 38, 42, 12.5, [(f"완전사례 주코호트 {r('Main analysis cohort').Total}", True),
                      (f"GPL14550 {r('Main analysis cohort').GPL14550} · GPL6480 {r('Main analysis cohort').GPL6480}", False),
                      ("PTX 표기 2명 포함 민감도 코호트 188 (126 · 62)", False)])
box(5, 6, 42, 9, [("최종 주분석 코호트 177", True), ("탐색 GPL14550 118 · 재현 GPL6480 59", False)], main=True)
box(56, 77, 40, 9, [("다른 진단군 362명 제외", True), ("간질성폐질환 254 · 대조군 108", False)], side=True)
box(56, 46.5, 40, 17.5, [("변수 불완전 34명 제외 (사유 중복)", True), ("DLCO 결측 23 · 폐기종 결측 9", False),
                     ("흡연 결측 3 · FEV1 결측 2", False), ("폐기종 문자 표기(PTX) 2", False)], side=True)
box(56, 18.5, 40, 13.5, [("배열 품질 이상치 9명 제외", True), ("GPL14550 7 · GPL6480 2", False),
                     ("배열 간 상관 중앙값 < Q1 − 3×IQR", False)], side=True)
kw = dict(arrowstyle="-|>", color=MUTED, lw=0.9, mutation_scale=9)
for y0, y1 in [(87.6, 75.4), (65.6, 50.9), (37.6, 15.4)]:
    ax.annotate("", xy=(26, y1), xytext=(26, y0), arrowprops=kw)
for y in [81.5, 55, 25]:
    ax.annotate("", xy=(55.6, y), xytext=(26, y), arrowprops=kw)
fig.savefig(out / "Figure1_cohort_flow.png"); plt.close(fig)

# ---------------- Figure 2: H1 (QC cohort) ----------------
h1 = pd.read_csv(A / "04_main/H1_primary_all_genes.csv")
fig, axs = plt.subplots(1, 3, figsize=(12, 3.9), gridspec_kw={"width_ratios": [1, 1, 1.05]})
for ax, pre, title, col in [(axs[0], "disc", "탐색 GPL14550 (n = 118)", BLUE), (axs[1], "rep", "재현 GPL6480 (n = 59)", ORANGE)]:
    x, y = h1[f"{pre}_beta_per10"], -np.log10(h1[f"{pre}_p"])
    sig = h1[f"{pre}_FDR"] < 0.05
    ax.scatter(x[~sig], y[~sig], s=3, c=GREY, lw=0, rasterized=True)
    ax.scatter(x[sig], y[sig], s=12, c=col, lw=0)
    ax.set_title(f"{title}: FDR < 0.05 유전자 {sig.sum()}개", loc="left", fontsize=9.5)
    ax.set_xlabel("DLCO +10%p당 log2 발현 변화 (보정)"); ax.set_ylabel("−log10 p")
    off = {"GPIHBP1": (4, 2), "LYRM1": (4, 2), "BTNL9": (4, 2), "BHLHE22": (-48, 8), "CYP1B1": (6, -4), "GDF15": (-38, -10)}
    for g, o in off.items():
        rr = h1[h1.symbol == g].iloc[0]
        xy = (rr[f"{pre}_beta_per10"], -np.log10(rr[f"{pre}_p"]))
        ax.scatter(*xy, s=14, facecolors="none", edgecolors=INK, lw=0.7)
        ax.annotate(g, xy, fontsize=7, color=INK, xytext=o, textcoords="offset points")
tag(axs[0], "A"); tag(axs[1], "B")
ax = axs[2]
top = h1.nsmallest(100, "disc_p")
ax.scatter(h1.disc_t, h1.rep_t, s=3, c=GREY, lw=0, rasterized=True, label="전체 11,394개")
ax.scatter(top.disc_t, top.rep_t, s=12, c=BLUE, lw=0, label="탐색 상위 100개")
meta = h1[h1.meta_FDR < 0.05]
ax.scatter(meta.disc_t, meta.rep_t, s=30, facecolors="none", edgecolors=INK, lw=0.8, label="메타 FDR < 0.05 (7개)")
moff = {"GPIHBP1": (5, -3), "NEB": (-28, 6), "NSMCE4A": (5, -9), "ZMAT3": (6, 0), "VWCE": (6, -10), "CYP1B1": (6, -3), "BHLHE22": (6, -4)}
for _, rr in meta.iterrows():
    ax.annotate(rr.symbol, (rr.disc_t, rr.rep_t), fontsize=7, xytext=moff.get(rr.symbol, (3, -8)), textcoords="offset points")
ax.axhline(0, color=MUTED, lw=0.6, ls=":"); ax.axvline(0, color=MUTED, lw=0.6, ls=":")
conc = (np.sign(top.disc_t) == np.sign(top.rep_t)).mean()
ax.set_title(f"플랫폼 간 t 통계량: 상위 100개 방향 일치 {conc*100:.0f}%", loc="left", fontsize=9.5)
ax.set_xlabel("탐색 t (GPL14550)"); ax.set_ylabel("재현 t (GPL6480)")
ax.legend(frameon=False, fontsize=7, loc="lower right"); tag(ax, "C")
fig.tight_layout(); fig.savefig(out / "Figure2_H1_volcano_concordance.png"); plt.close(fig)

# ---------------- Figure 3: replication diagnostics + QC ----------------
d1 = pd.read_csv(A / "03_diagnostics/diag1_split_half.csv")
d6 = pd.read_csv(A / "03_diagnostics/diag6_split_half_QC.csv")
lo_pre, hi_pre = [float(v.split("/")[0]) for v in d1.loc[1, ["q05", "q95"]]]
lo_q, hi_q = map(float, d6.loc[1, "conc_q05_q95"].split("-"))
rows = [("QC 전", "탐색 내부 split-half (64→61)", d1.loc[1, "conc_top500"], lo_pre, hi_pre, BLUE),
        ("QC 전", "교차 플랫폼 (125→61)", d1.loc[0, "conc_top500"], None, None, ORANGE),
        ("QC 후", "탐색 내부 split-half (59→59)", d6.loc[1, "conc_top500"], lo_q, hi_q, BLUE),
        ("QC 후", "교차 플랫폼 (118→59)", d6.loc[0, "conc_top500"], None, None, ORANGE)]
fig = plt.figure(figsize=(12, 4.2))
ax = fig.add_axes([0.07, 0.15, 0.36, 0.72])
ys = [3.6, 2.8, 1.2, 0.4]
for (grp, lab, m, lo, hi, c), y in zip(rows, ys):
    if lo is not None:
        ax.plot([lo, hi], [y, y], color=c, lw=2.2, solid_capstyle="round")
    ax.plot(m, y, "o" if lo is not None else "D", color=c, ms=7)
    ax.text(m, y + 0.22, f"{m:.2f}", ha="center", fontsize=8, color=INK)
ax.set_yticks(ys); ax.set_yticklabels([r[1] for r in rows])
ax.axhspan(2.3, 4.1, color=GRID, alpha=0.5, lw=0); ax.text(0.31, 3.95, "QC 전", fontsize=8.5, fontweight="bold", va="top")
ax.text(0.31, 1.75, "QC 후", fontsize=8.5, fontweight="bold", va="top")
ax.axvline(0.5, color=MUTED, lw=0.7, ls=":"); ax.set_xlim(0.3, 1.0); ax.set_ylim(0, 4.2)
ax.set_xlabel("상위 500개 유전자 효과 방향 일치율 (0.5 = 우연)")
ax.set_title("재현 일치율: split-half는 중앙값과 5–95% 범위 (200회)", loc="left", fontsize=9.5); tag(ax, "A")
ts = pd.read_csv(A / "05_raw_qc/tissue_content_scores_cohort186.csv")
cols = [("alveolar_AT2", "AT2 (폐포)"), ("alveolar_cap", "폐포 모세혈관"), ("cilia", "섬모"),
        ("airway_secretory", "기도 분비"), ("adipose", "지방"), ("lymphoid", "림프")]
ax = fig.add_axes([0.53, 0.15, 0.45, 0.72])
rng = np.random.default_rng(1)
for i, (c, lab) in enumerate(cols):
    a = ts.loc[~ts.qc_outlier, c]; b = ts.loc[ts.qc_outlier, c]
    ax.scatter(i + rng.uniform(-0.18, 0.18, len(a)), a, s=6, c=GREY, lw=0)
    ax.scatter(np.full(len(b), i) + rng.uniform(-0.1, 0.1, len(b)), b, s=22, c=ORANGE, edgecolors="white", lw=0.4)
ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[1] for c in cols])
ax.axhline(0, color=MUTED, lw=0.6, ls=":")
ax.set_ylabel("조직 구성 표지 점수 (z 평균)")
ax.scatter([], [], s=6, c=GREY, label="QC 통과 177"); ax.scatter([], [], s=22, c=ORANGE, label="QC 이상치 9")
ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=2)
ax.set_title("QC 이상치는 폐포 실질 표지가 낮고 비폐포 조직 표지가 높음", loc="left", fontsize=9.5); tag(ax, "B")
fig.savefig(out / "Figure3_replication_QC.png"); plt.close(fig)

# ---------------- Figure 4: H2 forest ----------------
h2 = pd.read_csv(A / "04_main/H2_genesets.csv")
mt = pd.read_csv(A / "04_main/H2_meta_exploratory.csv").set_index("set")
sets = [("H2_AT1", "AT1"), ("H2_gCap", "gCap"), ("H2_aCap", "aCap"), ("H2_AT2", "AT2")]
fig, ax = plt.subplots(figsize=(7.4, 4.6))
y = 0; yt, yl = [], []
for s, lab in sets:
    for p, c, nm in [("GPL14550", BLUE, "탐색"), ("GPL6480", ORANGE, "재현")]:
        rr = h2[(h2.platform == p) & (h2.set == s)].iloc[0]
        ax.plot([rr.beta_SD_per10_adj - 1.96 * rr.se, rr.beta_SD_per10_adj + 1.96 * rr.se], [y, y], color=c, lw=1.8)
        ax.plot(rr.beta_SD_per10_adj, y, "o", color=c, ms=5.5)
        ax.text(0.42, y, f"p = {rr.p_adj:.3f}", fontsize=7.5, va="center", color=MUTED)
        yt.append(y); yl.append(f"{lab} {nm}"); y -= 1
    m = mt.loc[s]
    lo, hi, b = m.meta_beta_per10 - 1.96 * m.meta_se, m.meta_beta_per10 + 1.96 * m.meta_se, m.meta_beta_per10
    ax.fill([lo, b, hi, b], [y, y + 0.28, y, y - 0.28], color=INK)
    ax.text(0.42, y, f"p = {m.meta_p:.4f} (Holm {m.meta_p_holm:.3f})", fontsize=7.5, va="center", color=INK)
    yt.append(y); yl.append(f"{lab} 메타"); y -= 1.6
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=8)
ax.axvline(0, color=MUTED, lw=0.7, ls=":"); ax.set_xlim(-0.2, 0.62)
ax.set_xlabel("DLCO +10%p당 GSVA 점수 변화 (SD), FEV1·폐기종·나이·성별·흡연 보정, 95% CI")
ax.set_title("H2: 사전 정의 폐포 세포유형 점수와 DLCO", loc="left", fontsize=10)
fig.savefig(out / "Figure4_H2_forest.png"); plt.close(fig)

# ---------------- Figure 5: H4 ----------------
h4 = pd.read_csv(A / "04_main/H4_prediction.csv")
labs = ["A: 임상변수", "B: A + H2 점수 4개", "C: A + 탐색 100유전자 서명"]
fig, ax = plt.subplots(figsize=(6.8, 2.8))
for i, (_, rr) in enumerate(h4.iterrows()):
    yy = 2 - i
    ax.plot([rr.R2_test_CI_low, rr.R2_test_CI_high], [yy, yy], color=ORANGE, lw=2)
    ax.plot(rr.R2_test, yy, "o", color=ORANGE, ms=7)
    ax.plot(rr.R2_train, yy, "|", color=BLUE, ms=14, mew=2)
    ax.text(0.80, yy, f"재현 R² {rr.R2_test:.3f}" + ("" if i == 0 else f"  (ΔR² {rr.dR2_test_vs_A:+.3f}, {rr.dR2_CI.replace(' to ', '~')})"),
            fontsize=7.5, va="center")
ax.set_yticks([2, 1, 0]); ax.set_yticklabels(labs)
ax.set_xlim(0.4, 0.8); ax.set_ylim(-0.6, 2.6); ax.set_xlabel("DLCO 설명력 R²")
ax.plot([], [], "|", color=BLUE, ms=10, mew=2, label="탐색 학습 R²")
ax.plot([], [], "o", color=ORANGE, label="재현 검증 R² (부트스트랩 95% CI)")
ax.legend(frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2)
fig.savefig(out / "Figure5_H4_prediction.png"); plt.close(fig)

# ---------------- Figure S5: array QC metrics ----------------
fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
for ax, p, c in [(axs[0], "GPL14550", BLUE), (axs[1], "GPL6480", ORANGE)]:
    q = pd.read_csv(Q / f"qc_array_{p}.csv")
    q1, q3 = q.medcor.quantile([.25, .75])
    ok = ~q.sample_id.isin(OUTL)
    ax.scatter(q.medcor[ok], q.PC_maha[ok], s=10, c=c, alpha=0.6, lw=0)
    ax.scatter(q.medcor[~ok], q.PC_maha[~ok], s=28, facecolors="none", edgecolors=INK, lw=1)
    for _, rr in q[~ok].iterrows():
        ax.annotate(rr.sample_id.replace("GSM", ""), (rr.medcor, rr.PC_maha), fontsize=6.5, xytext=(3, 2), textcoords="offset points")
    ax.axvline(q1 - 3 * (q3 - q1), color=MUTED, ls="--", lw=0.8)
    ax.set_title(f"{p} (n = {len(q)}): 이상치 {(~ok).sum()}개", loc="left", fontsize=9.5)
    ax.set_xlabel("다른 배열과의 상관 중앙값 (medcor)"); ax.set_ylabel("PC1–5 제곱 Mahalanobis 거리")
fig.tight_layout(); fig.savefig(out / "FigureS5_array_QC_metrics.png"); plt.close(fig)
print("figures written to", out)
