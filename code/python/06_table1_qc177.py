"""
GSE47460 COPD DLCO study - Step 6: Table 1 for the primary (array-QC) cohort, n = 177 (GPL14550 118 / GPL6480 59).

Usage: python 06_table1_qc177.py <01_cohort dir> <qc_array dir> <out dir>
  <qc_array dir> holds qc_array_GPL14550.csv / qc_array_GPL6480.csv (from 03b_array_qc_metrics.py),
  whose qc_outlier column defines the 9 excluded arrays.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

cdir, qdir, out = map(Path, sys.argv[1:4])
out.mkdir(parents=True, exist_ok=True)
d = pd.read_csv(cdir / "cohort_main_186.tsv", sep="\t")
qc = pd.concat([pd.read_csv(qdir / f"qc_array_{p}.csv") for p in ["GPL14550", "GPL6480"]])
d = d[~d.sample_id.isin(qc.sample_id[qc.qc_outlier])].copy()
d["sqrtEmph"] = np.sqrt(d.Emphysema)
assert len(d) == 177, len(d)
g1, g2 = d[d.platform == "GPL14550"], d[d.platform == "GPL6480"]

rows = [("N", len(d), len(g1), len(g2), "", "")]
for v, lab in [("Age", "나이, 세"), ("FEV1", "FEV1 %예측치 (기관지확장제 전)"), ("FVC", "FVC %예측치 (기관지확장제 전)"),
               ("DLCO", "DLCO %예측치"), ("Emphysema", "CT 폐기종 %LAA-950")]:
    f = lambda s: f"{s.mean():.1f} ± {s.std():.1f}"
    q = lambda s: f"{s.median():.1f} [{s.quantile(.25):.1f}–{s.quantile(.75):.1f}]"
    p = stats.mannwhitneyu(g1[v], g2[v]).pvalue
    rows.append((lab + ", 평균 ± SD", f(d[v]), f(g1[v]), f(g2[v]), f"{p:.3f}", "Mann-Whitney U"))
    rows.append(("   중앙값 [IQR]", q(d[v]), q(g1[v]), q(g2[v]), "", ""))
for v, lab, levels, names in [
        ("Sex", "성별", ["Male", "Female"], ["남성", "여성"]),
        ("Smoking", "흡연 상태", ["Current", "Former (ever >100 cig)", "Never"], ["현재 흡연", "과거 흡연", "비흡연"]),
        ("GOLD", "GOLD 병기", ["At Risk", "Mild", "Moderate", "Severe", "Very Severe", "Missing"],
         ["0 (At Risk)", "1 (Mild)", "2 (Moderate)", "3 (Severe)", "4 (Very Severe)", "결측"])]:
    ct = pd.crosstab(d[v].fillna("Missing"), d.platform).reindex(levels).dropna(how="all").fillna(0).astype(int)
    tab = ct.drop(index="Missing", errors="ignore")
    if tab.shape == (2, 2):
        p, test = stats.fisher_exact(tab.values).pvalue, "Fisher 정확검정"
    else:
        p, test = stats.chi2_contingency(tab.values)[1], "카이제곱 (기대빈도 5 미만 칸 있음)"
    rows.append((lab + ", n (%)", "", "", "", f"{p:.3f}", test))
    for lv, nm in zip(levels, names):
        if lv not in ct.index:
            continue
        a, b = ct.loc[lv, "GPL14550"], ct.loc[lv, "GPL6480"]
        rows.append((f"   {nm}", f"{a+b} ({(a+b)/len(d)*100:.1f})", f"{a} ({a/len(g1)*100:.1f})", f"{b} ({b/len(g2)*100:.1f})", "", ""))
t = pd.DataFrame(rows, columns=["변수", "전체", "GPL14550 (탐색)", "GPL6480 (재현)", "p", "검정"])
t.to_csv(out / "table1_qc177.csv", index=False, encoding="utf-8-sig")
print(t.to_string(index=False))
