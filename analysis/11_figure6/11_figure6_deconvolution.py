"""Step 11: Figure 6 (deconvolution) for the integrated manuscript.
(A) estimated cell proportions (%) ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current: beta per +10%p DLCO, discovery / replication / IVW meta.
(B) H2 GSVA scores before vs after adjustment for deconvolved non-alveolar composition (meta beta, 95% CI).
SEs are recovered from beta and two-sided p (t with n-7 df per platform; normal for meta).
Usage: python 11_figure6_deconvolution.py <10_deconvolution dir> <out dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D, out = map(Path, sys.argv[1:3]); out.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1f1f1e", "#6b6b68", "#e6e6e3"
plt.rcParams.update({"font.family": ["Malgun Gothic", "Noto Sans CJK KR", "Noto Sans CJK JP", "DejaVu Sans"], "axes.unicode_minus": False,
    "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "savefig.bbox": "tight"})
se_t = lambda b, p, df: abs(b) / stats.t.ppf(1 - p / 2, df)
se_z = lambda b, p: abs(b) / stats.norm.ppf(1 - p / 2)
pr = pd.read_csv(D / "deconv_proportions_vs_DLCO.csv").set_index("proportion")
h2 = pd.read_csv(D / "H5_deconv_H2.csv").set_index("set")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 4.6), gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.55))
rows = [("AlvCap", "폐포-모세혈관 합계"), ("AT1", "AT1"), ("AT2", "AT2"), ("CapA", "aCap (aerocyte)"), ("CapB", "gCap"),
        ("Macrophage", "대식세포"), ("Stromal", "기질세포"), ("Lymphoid", "림프구"), ("Airway", "기도상피")]
y = 0; yt, yl = [], []
for k, lab in rows:
    r = pr.loc[k]
    for b, p, df, c, dy in [(r.disc_beta_pct_per10, r.disc_p, 111, BLUE, 0.22), (r.rep_beta_pct_per10, r.rep_p, 52, ORANGE, 0)]:
        s = se_t(b, p, df); a1.plot([b - 1.96 * s, b + 1.96 * s], [y + dy, y + dy], color=c, lw=1.6); a1.plot(b, y + dy, "o", color=c, ms=4.5)
    m, s = r.meta_beta, se_z(r.meta_beta, r.meta_p); yy = y - 0.24
    a1.fill([m - 1.96 * s, m, m + 1.96 * s, m], [yy, yy + 0.15, yy, yy - 0.15], color=INK)
    a1.text(2.35, y, f"p = {r.meta_p:.4f}" if r.meta_p < 0.01 else f"p = {r.meta_p:.3f}", fontsize=7.5, va="center",
            color=INK if r.meta_p < 0.05 else MUTED, fontweight="bold" if r.meta_p < 0.05 else "normal")
    yt.append(y); yl.append(lab); y -= 1 if k != "CapB" else 1.5
a1.set_yticks(yt); a1.set_yticklabels(yl, fontsize=8.2); a1.axvline(0, color=MUTED, lw=0.7, ls=":"); a1.set_xlim(-2.6, 2.3)
a1.text(2.35, 0.75, "메타 p", fontsize=7.5, color=MUTED)
a1.set_xlabel("DLCO +10%p당 추정 세포 비율 변화 (%p), 95% CI")
a1.set_title("추정 세포 구성비와 잔여 DLCO", loc="left", fontsize=10)
a1.plot([], [], "o-", color=BLUE, label="탐색 (n = 118)"); a1.plot([], [], "o-", color=ORANGE, label="재현 (n = 59)")
a1.plot([], [], "D", color=INK, label="역분산 가중 메타"); a1.legend(frameon=False, fontsize=7.5, loc="lower left")
a1.text(-0.12, 1.04, "A", transform=a1.transAxes, fontsize=12, fontweight="bold")
y = 0; yt, yl = [], []
for k in ["AT1", "gCap", "aCap", "AT2"]:
    r = h2.loc[k]
    for b, p, c, lab, dy in [(r.meta_beta, r.meta_p, GRID, "보정 전", 0.18), (r.meta_beta_adj, r.meta_p_adj, INK, "구성비 보정 후", -0.18)]:
        s = se_z(b, p); a2.plot([b - 1.96 * s, b + 1.96 * s], [y + dy, y + dy], color="#9a9a96" if c == GRID else INK, lw=1.8)
        a2.plot(b, y + dy, "o", color="#9a9a96" if c == GRID else INK, ms=5, mfc="white" if c == GRID else INK)
        a2.text(0.33, y + dy, f"p = {p:.4f}" if p < 0.01 else f"p = {p:.3f}", fontsize=7.5, va="center", color=MUTED)
    yt.append(y); yl.append(k); y -= 1
a2.set_yticks(yt); a2.set_yticklabels(yl, fontsize=8.5); a2.axvline(0, color=MUTED, lw=0.7, ls=":"); a2.set_xlim(-0.12, 0.32)
a2.plot([], [], "o-", color="#9a9a96", mfc="white", label="보정 전 (주분석)"); a2.plot([], [], "o-", color=INK, label="비폐포 구성비 보정 후")
a2.legend(frameon=False, fontsize=7.5, loc="lower right")
a2.set_xlabel("DLCO +10%p당 GSVA 점수 변화 (SD), 메타 95% CI")
a2.set_title("H2 점수의 구성비 보정 전후", loc="left", fontsize=10)
a2.text(-0.12, 1.04, "B", transform=a2.transAxes, fontsize=12, fontweight="bold")
fig.savefig(out / "Figure6_deconvolution.png"); print("ok")
