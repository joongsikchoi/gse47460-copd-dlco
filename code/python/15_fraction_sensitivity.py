"""Sensitivity analyses of the alveolar-capillary fraction (added for submission to Respiratory Medicine).
Usage: python 15_fraction_sensitivity.py <analysis_dir> <output_dir>
Inputs: <analysis_dir>/01_cohort/cohort_qc_177.tsv, <analysis_dir>/10_deconvolution/deconv_proportions_186.csv
Output: <output_dir>/fraction_sensitivity.csv (Supplementary Table S16)
Model per platform: y ~ DLCO/10 + FEV1 + sqrt(emphysema) + age + male + current smoking; inverse-variance fixed-effect meta.
"""
import sys, pandas as pd, numpy as np, statsmodels.formula.api as smf
from scipy import stats
A = sys.argv[1] if len(sys.argv) > 1 else "analysis"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{A}/15_fraction_sensitivity"
import os; os.makedirs(OUT, exist_ok=True)
c = pd.read_csv(f"{A}/01_cohort/cohort_qc_177.tsv", sep="\t")
d = pd.read_csv(f"{A}/10_deconvolution/deconv_proportions_186.csv")
df = c.merge(d, on=["sample_id", "platform"])
df["Male"] = (df.Sex == "Male").astype(int); df["Current"] = (df.Smoking == "Current").astype(int)
df["sq"] = np.sqrt(df.Emphysema); df["D10"] = df.DLCO / 10
df["upper"] = df.title.str.extract(r"LT\d{6}(\w\w)")[0].isin(["RU", "LU"]).astype(int)
alv = df[["AT1", "AT2", "CapA", "CapB"]].sum(axis=1)
imm = df[["Monocyte", "DC", "Mast", "T_NK", "B_Plasma", "Macrophage"]].sum(axis=1)
stro = df[["Fibroblast", "SMC_Pericyte", "VE_other", "Lymphatic"]].sum(axis=1)
air = df[["Airway_epi", "Ciliated"]].sum(axis=1)
df["alv_pct"] = alv * 100
df["alv_noAT2_pct"] = df[["AT1", "CapA", "CapB"]].sum(axis=1) * 100
df["logit_alv"] = np.log(alv / (1 - alv))
df["lr_imm"] = np.log(alv / imm); df["lr_str"] = np.log(alv / stro); df["lr_air"] = np.log(alv / air)
F = "D10+FEV1+sq+Age+Male+Current"
gold0 = df.GOLD.astype(str).str.contains("0|at risk", case=False)
rows = []
def meta(label, y, data=df, extra=""):
    b, s = [], []
    for p in ["GPL14550", "GPL6480"]:
        m = smf.ols(f"{y}~{F}{extra}", data[data.platform == p]).fit(); b.append(m.params.D10); s.append(m.bse.D10)
    b, s = np.array(b), np.array(s); w = 1 / s**2
    mb = (w * b).sum() / w.sum(); ms = 1 / np.sqrt(w.sum()); p = 2 * stats.norm.sf(abs(mb / ms))
    Q = (w * (b - mb) ** 2).sum(); Qp = stats.chi2.sf(Q, 1)
    rows.append(dict(analysis=label, n=f"{(data.platform=='GPL14550').sum()}/{(data.platform=='GPL6480').sum()}",
                     disc=f"{b[0]:.3f}", rep=f"{b[1]:.3f}", meta=f"{mb:.3f} [{mb-1.96*ms:.3f}, {mb+1.96*ms:.3f}]",
                     p=f"{p:.2g}", Q_p=f"{Qp:.2f}"))
    return mb, s
mb, s = meta("Main: alveolar-capillary fraction (AT1+AT2+aCap+gCap), %", "alv_pct")
meta("Without AT2 (AT1+aCap+gCap), %", "alv_noAT2_pct")
meta("GOLD 0 excluded, %", "alv_pct", data=df[~gold0])
meta("Adjusted for upper-lobe sample, %", "alv_pct", extra="+upper")
meta("Logit of the fraction", "logit_alv")
meta("Log-ratio to immune cells", "lr_imm")
meta("Log-ratio to stromal and vascular cells", "lr_str")
meta("Log-ratio to airway epithelial cells", "lr_air")
out = pd.DataFrame(rows); print(out.to_string(index=False))
m = smf.ols(f"alv_pct~{F}", df[df.platform == "GPL6480"]).fit(); lo, hi = m.conf_int().loc["D10"]
print(f"\nReplication-platform estimate {m.params.D10:.3f}, 95% CI {lo:.2f} to {hi:.2f} (discovery estimate inside: {lo < 1.016 < hi})")
out.to_csv(f"{OUT}/fraction_sensitivity.csv", index=False)
