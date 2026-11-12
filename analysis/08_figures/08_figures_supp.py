"""Step 8: final figures and supplementary table for the QC-177 manuscript.
Figure 1: cohort flow (582 -> 177, by platform). Figure 3: split-half vs cross-platform concordance, before/after array QC.
Figure 4: H2 forest with per-platform and inverse-variance meta estimates. Table S1: QC-excluded arrays.
Usage: python 08_figures_supp.py <gse dir with out/, out1q/, out3/, out4/, out5/, out8/> <out dir>"""
import sys, os, numpy as np, pandas as pd, statsmodels.formula.api as smf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
g, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#222222", "#6b6b6b", "#e6e6e3"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, name):
    for ext in ("png", "pdf"): fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---------------- Figure 1: cohort flow ----------------
flow = pd.read_csv(os.path.join(g, "out1q/cohort_flow_qc_177.csv")).set_index("Step")
miss = flow.loc[[i for i in flow.index if "Missing" in i or "non-numeric" in i]]
boxes = [("All GSE47460 arrays", 582, 429, 153), ("COPD", 220, 145, 75), ("Complete cases", 188, 126, 62),
         ("Numeric emphysema", 186, 125, 61), ("Main analysis (after array QC)", 177, 118, 59)]
side = ["Non-COPD (ILD, controls): 362 (284 / 78)",
        "Missing ≥1 variable: 32 (19 / 13)\nDLCO 23, emphysema 9, smoking 3, FEV1 2 (overlapping)",
        "Emphysema value flagged \"(PTX)\": 2 (1 / 1)",
        "Array QC outliers*: 9 (7 / 2)"]
assert boxes[1][1] - boxes[2][1] == 32 and boxes[0][1] - boxes[1][1] == 362
fig, ax = plt.subplots(figsize=(7.2, 6.2)); ax.set_xlim(0, 11.2); ax.set_ylim(0, 10.4); ax.axis("off")
ys = [9.4, 7.3, 5.2, 3.1, 1.0]
for (t, n, a, b), y in zip(boxes, ys):
    last = t.startswith("Main")
    ax.add_patch(FancyBboxPatch((0.3, y - 0.55), 4.6, 1.1, boxstyle="round,pad=0.02,rounding_size=0.12",
                 fc="#eef4fc" if last else "white", ec=BLUE if last else MUTED, lw=1.4 if last else 0.9))
    ax.text(2.6, y + 0.18, f"{t}\nn = {n}", ha="center", va="center", fontsize=9, color=INK, weight="bold" if last else "normal")
    ax.text(2.6, y - 0.36, f"GPL14550 {a}  |  GPL6480 {b}", ha="center", va="center", fontsize=8, color=MUTED)
for i, s in enumerate(side):
    y0, y1 = ys[i] - 0.55, ys[i + 1] + 0.55; ym = (y0 + y1) / 2
    ax.annotate("", xy=(2.6, y1), xytext=(2.6, y0), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9))
    ax.plot([2.6, 5.3], [ym, ym], color=MUTED, lw=0.9)
    ax.add_patch(FancyBboxPatch((5.3, ym - 0.42), 5.85, 0.84, boxstyle="round,pad=0.02,rounding_size=0.1", fc="#f6f6f4", ec=GRID, lw=0.8))
    ax.text(5.45, ym, "Excluded — " + s, ha="left", va="center", fontsize=7.3, color=INK)
ax.text(0.3, -0.2, "Counts are total (GPL14550 / GPL6480). *Median correlation with other arrays below Q1 − 3×IQR within platform;\n"
        "8 of 9 had low alveolar or high non-alveolar tissue marker scores (Table S1). GPL14550 = discovery, GPL6480 = replication; platforms not merged.",
        fontsize=6.8, color=MUTED, va="top")
save(fig, "Figure1_cohort_flow")

# ---------------- Figure 3: split-half vs cross-platform ----------------
d = pd.read_csv(os.path.join(g, "out8/split_half_draws.csv")); x = pd.read_csv(os.path.join(g, "out8/cross_platform_top500.csv")).set_index("cohort")
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
s3 = []
for ax, (c, lab) in zip(axs, [("preQC_186", "A  Before array QC (n = 125 / 61)"), ("QC_177", "B  After array QC (n = 118 / 59)")]):
    v = d.loc[d.cohort == c, "conc"] * 100; cp = x.loc[c, "conc"] * 100; frac = (v <= cp).mean()
    ax.hist(v, bins=np.arange(40, 100, 2.5), color=BLUE, ec="white", lw=1)
    ax.axvline(cp, color=ORANGE, lw=2); ax.axvline(50, color=MUTED, lw=0.8, ls=(0, (3, 3)))
    ax.set_ylim(0, 37)
    pre = c == "preQC_186"
    ax.text(cp + 1.2 if pre else cp - 1.2, 35, f"Cross-platform\n{cp:.0f}%", ha="left" if pre else "right", va="top",
            fontsize=7.5, color=INK, bbox=dict(fc="white", ec="none", pad=1.5))
    ax.text(97, 35, f"Split-half median {v.median():.0f}%\n(5–95%: {v.quantile(.05):.0f}–{v.quantile(.95):.0f}%)",
            ha="right", va="top", fontsize=7.5, color=INK)
    ax.set_title(lab, loc="left", fontsize=8.5, color=INK); ax.set_xlabel("Direction concordance of top-500 genes (%)"); ax.set_xlim(40, 100)
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    s3.append(dict(cohort=c, cross_platform=round(cp, 1), split_median=round(v.median(), 1), split_q05=round(v.quantile(.05), 1),
                   split_q95=round(v.quantile(.95), 1), frac_splits_at_or_below_cross=round(frac, 3)))
axs[0].set_ylabel("Random split-halves (of 200)")
fig.text(0.01, -0.06, "Blue: discovery cohort split at random into a training half and a test half the size of the replication cohort. "
         "Orange: observed GPL14550 → GPL6480 concordance. Dashed: chance (50%).", fontsize=6.8, color=MUTED)
fig.tight_layout(); save(fig, "Figure3_split_half_vs_cross_platform")
pd.DataFrame(s3).to_csv(os.path.join(out, "Figure3_summary.csv"), index=False)

# ---------------- Figure 4: H2 forest ----------------
h2 = pd.read_csv(os.path.join(g, "out4/H2_genesets.csv")); m = pd.read_csv(os.path.join(g, "out4/H2_meta_exploratory.csv"))
order = ["H2_aCap", "H2_gCap", "H2_AT1", "H2_AT2"]; names = {"H2_aCap": "Aerocyte capillary (aCap)", "H2_gCap": "General capillary (gCap)",
                                                          "H2_AT1": "Alveolar type 1 (AT1)", "H2_AT2": "Alveolar type 2 (AT2)"}
fig, ax = plt.subplots(figsize=(7.2, 4.4)); y = 0; ticks = []; rows = []
for s in order[::-1]:
    mm = m[m.set == s].iloc[0]
    for lab, b, se, p, col, mk in [("Combined (IVW)", mm.meta_beta_per10, mm.meta_se, mm.meta_p, INK, "D"),
                                   ("GPL6480 (replication)",) + tuple(h2[(h2.set == s) & (h2.platform == "GPL6480")][["beta_SD_per10_adj", "se", "p_adj"]].iloc[0]) + (ORANGE, "o"),
                                   ("GPL14550 (discovery)",) + tuple(h2[(h2.set == s) & (h2.platform == "GPL14550")][["beta_SD_per10_adj", "se", "p_adj"]].iloc[0]) + (BLUE, "o")]:
        ax.plot([b - 1.96 * se, b + 1.96 * se], [y, y], color=col, lw=2, solid_capstyle="round")
        ax.plot(b, y, mk, color=col, ms=6.5 if mk == "D" else 6, mec="white", mew=1.2)
        ax.text(0.52, y, f"{b:.2f} ({b - 1.96*se:.2f}, {b + 1.96*se:.2f})", va="center", fontsize=7.5, color=INK)
        ax.text(0.80, y, f"{p:.3f}" if p >= 0.001 else "<0.001", va="center", fontsize=7.5, color=INK)
        ticks.append((y, lab)); rows.append(dict(set=names[s], estimate=lab, beta_SD_per10pct=round(b, 3), ci_low=round(b - 1.96*se, 3), ci_high=round(b + 1.96*se, 3), p=p))
        y += 1
    ax.text(-0.24, y - 0.2, names[s], fontsize=8.5, color=INK, weight="bold", va="bottom"); y += 1.1
ax.axvline(0, color=MUTED, lw=0.8)
ax.set_yticks([t[0] for t in ticks]); ax.set_yticklabels([t[1] for t in ticks], fontsize=7.5)
ax.set_xlim(-0.24, 0.9); ax.set_xticks([-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4]); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
ax.text(0.52, y - 0.4, "β (95% CI)", fontsize=7.5, weight="bold", color=INK); ax.text(0.80, y - 0.4, "p", fontsize=7.5, weight="bold", color=INK)
ax.set_xlabel("SD change in gene-set score per +10 %-points DLCO\n(adjusted for FEV1, √emphysema, age, sex, current smoking)", fontsize=8)
ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.spines["bottom"].set_bounds(-0.2, 0.4)
fig.text(0.01, -0.09, "GSVA scores of MSigDB C8 Travaglini lung cell-type sets, z-scored within platform; QC cohort 118 / 59.\n"
         "Combined = fixed-effect inverse-variance meta-analysis of the two platform estimates (exploratory; data not merged).", fontsize=6.8, color=MUTED)
save(fig, "Figure4_H2_forest_meta"); pd.DataFrame(rows).to_csv(os.path.join(out, "Figure4_data.csv"), index=False)

# ---------------- Table S1: QC-excluded arrays ----------------
c = pd.read_csv(os.path.join(g, "out/cohort_main_186.tsv"), sep="\t"); c["sqrtEmph"] = np.sqrt(c.Emphysema)
c["Male"] = (c.Sex == "Male").astype(int); c["Current"] = (c.Smoking == "Current").astype(int)
c["resid_DLCO"] = np.nan
for p, s in c.groupby("platform"):
    c.loc[s.index, "resid_DLCO"] = smf.ols("DLCO ~ FEV1 + sqrtEmph + Age + Male + Current", s).fit().resid
qa = pd.concat([pd.read_csv(os.path.join(g, f"out3/qc_array_{p}.csv")) for p in ["GPL14550", "GPL6480"]])
t = pd.read_csv(os.path.join(g, "out5/tissue_content_scores_cohort186.csv"))
t = t.merge(qa[["sample_id", "medcor", "z_medcor"]], on="sample_id").merge(c[["sample_id", "title", "GOLD", "DLCO", "resid_DLCO"]], on="sample_id")
cols = ["alveolar_AT2", "alveolar_cap", "airway_secretory", "cilia", "adipose", "mesothelium", "lymphoid", "blood_RBC"]
ex = t[t.qc_outlier].copy()
def note(r):
    n = []
    if r.alveolar_AT2 < -3: n.append("very low alveolar")
    if r.adipose > 2: n.append("adipose")
    if r.mesothelium > 1.5: n.append("mesothelium")
    if r.lymphoid > 1.8: n.append("lymphoid")
    if r.airway_secretory > 1 and r.cilia > 1.5: n.append("airway")
    if r.blood_RBC < -2.5: n.append("low blood")
    if str(r.title).endswith("_ILD"): n.append("GEO title says ILD")
    return "; ".join(n) if n else "low inter-array correlation only"
ex["pattern"] = ex.apply(note, axis=1)
ref = t[~t.qc_outlier]
refrow = {"sample_id": "Retained 177, median [IQR]", "platform": "", "GOLD": "", "DLCO": f"{ref.DLCO.median():.0f}",
          "resid_DLCO": f"{ref.resid_DLCO.median():.1f} [{ref.resid_DLCO.quantile(.25):.1f}, {ref.resid_DLCO.quantile(.75):.1f}]",
          "medcor": f"{ref.medcor.median():.3f}", "z_medcor": "", "pattern": ""}
for k in cols: refrow[k] = f"{ref[k].median():.2f} [{ref[k].quantile(.25):.2f}, {ref[k].quantile(.75):.2f}]"
keep = ["sample_id", "platform", "GOLD", "DLCO", "resid_DLCO", "medcor", "z_medcor"] + cols + ["pattern"]
ex = ex.sort_values(["platform", "alveolar_AT2"], ascending=[False, True])[keep]
ex["resid_DLCO"] = ex.resid_DLCO.round(1)
tab = pd.concat([ex.astype(object), pd.DataFrame([refrow])[keep]], ignore_index=True)
tab.to_csv(os.path.join(out, "TableS1_QC_excluded_arrays.csv"), index=False)
print(tab.to_string()); print(pd.DataFrame(s3).to_string())
print("excluded GPL14550 mean resid DLCO:", ex[ex.platform == "GPL14550"].resid_DLCO.mean().round(1))
