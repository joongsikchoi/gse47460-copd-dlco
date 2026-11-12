"""
# Box panels: Mann-Whitney U p-values (platform comparison).
GSE47460 COPD DLCO study - Step 1: analysis cohort, Table 1, clinical-variable distributions by platform.

Main cohort (decision 2026-10-20): COPD with complete DLCO %pred, pre-BD FEV1 %pred, CT %emphysema (LAA-950),
age, sex, smoking status. Emphysema values that are non-numeric (e.g. "71.6 (PTX)", pneumothorax-flagged) are
excluded from the main analysis (n=186; GPL14550 125 / GPL6480 61) and kept for a sensitivity cohort (n=188).

Usage: python 01_cohort_table1.py <path to prepared/all_sample_metadata.tsv> <output dir>
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

meta_path, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)

d = pd.read_csv(meta_path, sep="\t")
c = d[d["disease state"].str.startswith("Chronic Obstructive")].copy()

ren = {"%predicted dlco": "DLCO", "%predicted fev1 (pre-bd)": "FEV1", "%predicted fvc (pre-bd)": "FVC",
       "%emphysema (f-950)": "emph_raw", "age": "Age", "Sex": "Sex", "smoker?": "Smoking", "gold stage": "GOLD"}
c = c.rename(columns=ren)
for v in ["DLCO", "FEV1", "FVC", "Age"]:
    c[v] = pd.to_numeric(c[v], errors="coerce")
c["Emphysema"] = pd.to_numeric(c["emph_raw"], errors="coerce")
c["emph_flag"] = c["emph_raw"].notna() & c["Emphysema"].isna()          # e.g. "71.6 (PTX)"
c["Emphysema_lenient"] = c["Emphysema"].fillna(
    c["emph_raw"].astype(str).str.extract(r"^\s*([\d.]+)")[0].astype(float))

c["Sex"] = c["Sex"].str.replace(r"^\d-", "", regex=True)
c["Smoking"] = c["Smoking"].map({"1-Current": "Current", "2-Ever (>100)": "Former (ever >100 cig)", "3-Never": "Never"})
c["GOLD"] = c["GOLD"].str.replace(r"^\d-", "", regex=True).str.replace(" COPD", "")

base = ["DLCO", "FEV1", "Age", "Sex", "Smoking"]
main = c.dropna(subset=base + ["Emphysema"]).copy()
sens = c.dropna(subset=base + ["Emphysema_lenient"]).copy()
main["cohort"] = "main_186"
sens_only = sens[~sens.sample_id.isin(main.sample_id)]

# ---- cohort file + flow ----
cols = ["sample_id", "platform", "title", "Age", "Sex", "Smoking", "GOLD", "FEV1", "FVC", "DLCO", "Emphysema", "emph_raw"]
main[cols].to_csv(out / "cohort_main_186.tsv", sep="\t", index=False)
sens[cols].assign(Emphysema=sens["Emphysema_lenient"]).to_csv(out / "cohort_sensitivity_188.tsv", sep="\t", index=False)

flow = [("All GSE47460 samples", len(d), (d.platform == "GPL14550").sum(), (d.platform == "GPL6480").sum()),
        ("COPD", len(c), (c.platform == "GPL14550").sum(), (c.platform == "GPL6480").sum())]
for lab, sub in [("Missing DLCO", c.DLCO.isna()), ("Missing FEV1 (pre-BD)", c.FEV1.isna()),
                 ("Missing emphysema", c.emph_raw.isna()), ("Missing smoking", c.Smoking.isna()),
                 ("Emphysema non-numeric (PTX flag)", c.emph_flag)]:
    flow.append(("  - " + lab + " (overlapping)", sub.sum(), (sub & (c.platform == "GPL14550")).sum(),
                 (sub & (c.platform == "GPL6480")).sum()))
flow += [("Main analysis cohort", len(main), (main.platform == "GPL14550").sum(), (main.platform == "GPL6480").sum()),
         ("Sensitivity cohort (+PTX-flagged)", len(sens), (sens.platform == "GPL14550").sum(), (sens.platform == "GPL6480").sum())]
pd.DataFrame(flow, columns=["Step", "Total", "GPL14550", "GPL6480"]).to_csv(out / "cohort_flow.csv", index=False)


# ---- Table 1 ----
def table1(df):
    g1, g2 = df[df.platform == "GPL14550"], df[df.platform == "GPL6480"]
    rows = [("N", f"{len(df)}", f"{len(g1)}", f"{len(g2)}", "", "")]
    for v, lab in [("Age", "Age, years"), ("FEV1", "FEV1 %pred (pre-BD)"), ("FVC", "FVC %pred (pre-BD)"),
                   ("DLCO", "DLCO %pred"), ("Emphysema", "CT emphysema %LAA-950")]:
        f = lambda s: f"{s.mean():.1f} ± {s.std():.1f}"
        q = lambda s: f"{s.median():.1f} [{s.quantile(.25):.1f}, {s.quantile(.75):.1f}]"
        a, b = g1[v].dropna(), g2[v].dropna()
        p = stats.mannwhitneyu(a, b).pvalue
        miss = df[v].isna().sum()
        rows.append((lab + (f" (missing {miss})" if miss else ""), f(df[v].dropna()), f(a), f(b), f"{p:.3f}", "Mann-Whitney U"))
        rows.append(("   median [IQR]", q(df[v].dropna()), q(a), q(b), "", ""))
    for v, lab, levels in [("Sex", "Sex", ["Male", "Female"]),
                           ("Smoking", "Smoking", ["Current", "Former (ever >100 cig)", "Never"]),
                           ("GOLD", "GOLD stage", ["At Risk", "Mild", "Moderate", "Severe", "Very Severe"])]:
        ct = pd.crosstab(df[v].fillna("Missing"), df.platform)
        ct = ct.reindex([l for l in levels + ["Missing"] if l in ct.index])
        tab = ct.drop(index="Missing", errors="ignore")
        p = stats.fisher_exact(tab.values).pvalue if tab.shape == (2, 2) else stats.chi2_contingency(tab.values)[1]
        test = "Fisher exact" if tab.shape == (2, 2) else "Chi-square (sparse cells; interpret with caution)"
        rows.append((lab + ", n (%)", "", "", "", f"{p:.3f}", test))
        for l in ct.index:
            n1, n2 = ct.loc[l].get("GPL14550", 0), ct.loc[l].get("GPL6480", 0)
            rows.append((f"   {l}", f"{n1+n2} ({(n1+n2)/len(df)*100:.1f})", f"{n1} ({n1/len(g1)*100:.1f})",
                         f"{n2} ({n2/len(g2)*100:.1f})", "", ""))
    return pd.DataFrame(rows, columns=["Variable", "Overall", "GPL14550 (discovery)", "GPL6480 (replication)", "p", "Test"])

t_main = table1(main)
t_sens = table1(sens.assign(Emphysema=sens["Emphysema_lenient"]))
t_main.to_csv(out / "table1_main_186.csv", index=False, encoding="utf-8-sig")
t_sens.to_csv(out / "table1_sensitivity_188.csv", index=False, encoding="utf-8-sig")

# ---- correlations among key covariates (multicollinearity preview) ----
corr = {p: main[main.platform == p][["DLCO", "FEV1", "Emphysema", "Age"]].corr(method="spearman").round(2)
        for p in ["GPL14550", "GPL6480"]}
with open(out / "covariate_spearman.txt", "w") as fh:
    for p, m in corr.items():
        fh.write(f"## {p}\n{m.to_string()}\n\n")

# ---- Figure ----
COL = {"GPL14550": "#2a78d6", "GPL6480": "#eb6834"}
INK, MUTED, GRID = "#1f1f1e", "#6b6b68", "#e6e6e3"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
plats = ["GPL14550", "GPL6480"]
labs = {p: f"{p}\n(n={(main.platform == p).sum()})" for p in plats}
fig, axes = plt.subplots(2, 4, figsize=(13, 6.6))
rng = np.random.default_rng(1)
cont = [("DLCO", "DLCO %pred"), ("FEV1", "FEV1 %pred"), ("Emphysema", "Emphysema %LAA-950"), ("Age", "Age, years")]
for ax, (v, lab) in zip(axes[0], cont):
    for i, p in enumerate(plats):
        y = main.loc[main.platform == p, v].values
        ax.boxplot(y, positions=[i], widths=0.5, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="none", edgecolor=MUTED, linewidth=1),
                   medianprops=dict(color=INK, linewidth=2), whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(y)), y, s=12, color=COL[p], alpha=0.65, linewidths=0)
    pval = stats.mannwhitneyu(*[main.loc[main.platform == p, v] for p in plats]).pvalue
    ax.set_xticks([0, 1], [labs[p] for p in plats]); ax.set_title(lab, color=INK, fontsize=10, loc="left")
    ax.set_title(f"p={pval:.2f}", color=MUTED, fontsize=8, loc="right")
    ax.yaxis.grid(True, color=GRID, linewidth=0.8); ax.set_axisbelow(True)

cat = [("Sex", "Sex", ["Male", "Female"]), ("Smoking", "Smoking", ["Current", "Former (ever >100 cig)", "Never"]),
       ("GOLD", "GOLD stage", ["At Risk", "Mild", "Moderate", "Severe", "Very Severe", "Missing"])]
for ax, (v, lab, levels) in zip(axes[1][:3], cat):
    ct = pd.crosstab(main[v].fillna("Missing"), main.platform, normalize="columns").reindex(levels).fillna(0) * 100
    x = np.arange(len(levels)); w = 0.38
    for j, p in enumerate(plats):
        ax.bar(x + (j - 0.5) * (w + 0.02), ct[p], width=w, color=COL[p], label=p)
    short = [l.replace("Former (ever >100 cig)", "Former") for l in levels]
    ax.set_xticks(x, short, rotation=0 if len(levels) < 4 else 30, ha="center" if len(levels) < 4 else "right")
    ax.set_ylabel("% within platform"); ax.set_title(lab, color=INK, fontsize=10, loc="left")
    ax.yaxis.grid(True, color=GRID, linewidth=0.8); ax.set_axisbelow(True)
axes[1][0].legend(frameon=False, fontsize=8)

ax = axes[1][3]
for p in plats:
    s = main[main.platform == p]
    ax.scatter(s.FEV1, s.DLCO, s=14, color=COL[p], alpha=0.7, linewidths=0, label=p)
ax.set_xlabel("FEV1 %pred"); ax.set_ylabel("DLCO %pred"); ax.set_title("DLCO vs FEV1", color=INK, fontsize=10, loc="left")
r = stats.spearmanr(main.FEV1, main.DLCO)[0]
ax.text(0.02, 0.98, f"Spearman ρ={r:.2f} (pooled)", transform=ax.transAxes, va="top", color=MUTED, fontsize=8)
ax.legend(frameon=False, fontsize=8, loc="lower right"); ax.grid(True, color=GRID, linewidth=0.8); ax.set_axisbelow(True)

fig.suptitle(f"GSE47460 COPD analysis cohort (n={len(main)}): clinical variables by array platform",
             x=0.01, ha="left", fontsize=12, color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(out / "fig1_clinical_by_platform.png", dpi=300)
fig.savefig(out / "fig1_clinical_by_platform.pdf")
print(t_main.to_string(index=False)); print(open(out / "covariate_spearman.txt").read())
print(pd.read_csv(out / "cohort_flow.csv").to_string(index=False))
