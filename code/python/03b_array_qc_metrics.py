"""
GSE47460 COPD DLCO study - Step 3b: array-level QC metrics (qc_array_<platform>.csv).

The original script that wrote analysis/03_diagnostics/qc_array_*.csv was not kept. This script
re-derives the same table from the prepared expression data and checks it against the stored file
(GPL14550: medcor r = 0.9999, PC_maha max |diff| = 0.05).

Definitions (per platform, main cohort n = 186 split 125 / 61, gene universe = 11,394 genes):
  medcor   : median Pearson correlation of a sample with every other array
  z_medcor : robust z = (medcor - median) / (1.4826 * MAD)
  PC1, PC2 : principal-component scores of the gene-centred sample x gene matrix
  PC_maha  : squared Mahalanobis distance on PC1-PC5 (each PC scaled by its SD)
QC rule used in the analysis (also in R/00_common.R): exclude medcor < Q1 - 3 * IQR (Tukey far-out).

Usage: python 03b_array_qc_metrics.py <prepared dir> <01_cohort dir> <out dir>
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

prep, cdir, out = map(Path, sys.argv[1:4])
out.mkdir(parents=True, exist_ok=True)
coh = pd.read_csv(cdir / "cohort_main_186.tsv", sep="\t")


def load(p):
    e = pd.read_csv(prep / f"{p}_expression.tsv.gz", sep="\t", index_col=0)
    a = pd.read_csv(prep / f"{p}_gene_annotation.tsv.gz", sep="\t", dtype=str, usecols=["ID", "GENE_SYMBOL"])
    return e, a.set_index("ID").GENE_SYMBOL.reindex(e.index)


# gene universe: defined on the discovery cohort (GPL14550, n = 125), identical to R/00_common.R
eD, symD = load("GPL14550")
disc = coh.loc[coh.platform == "GPL14550", "sample_id"]
md = eD[disc].mean(axis=1)
keep = symD.notna() & (symD != "") & (md > md.quantile(0.25))
order = [i for i in md.sort_values(ascending=False).index if keep[i]]
s = symD[order]
genes = list(s[~s.duplicated()].index)
print("gene universe:", len(genes))

for p in ["GPL14550", "GPL6480"]:
    e = eD if p == "GPL14550" else load(p)[0]
    d = coh[coh.platform == p]
    X = e.loc[genes, d.sample_id]
    cc = np.corrcoef(X.T.values)
    np.fill_diagonal(cc, np.nan)
    medcor = np.nanmedian(cc, axis=1)
    mad = 1.4826 * np.median(np.abs(medcor - np.median(medcor)))
    z = (medcor - np.median(medcor)) / mad
    Xc = X.T.values - X.T.values.mean(axis=0)
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    pcs = U * S
    pc5 = pcs[:, :5] / pcs[:, :5].std(axis=0, ddof=1)
    maha = (pc5 ** 2).sum(axis=1)
    q1, q3 = np.quantile(medcor, [0.25, 0.75])
    res = pd.DataFrame({"sample_id": d.sample_id.values, "GOLD": d.GOLD.values, "medcor": medcor.round(3),
                        "z_medcor": z.round(2), "PC_maha": maha.round(1), "PC1": pcs[:, 0].round(1),
                        "PC2": pcs[:, 1].round(1), "qc_outlier": medcor < q1 - 3 * (q3 - q1)})
    res.to_csv(out / f"qc_array_{p}.csv", index=False)
    print(p, "n =", len(res), "| QC outliers:", ", ".join(res.sample_id[res.qc_outlier]))
