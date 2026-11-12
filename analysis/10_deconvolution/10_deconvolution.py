"""Step 10: reference-based cell-type deconvolution (H5 with a single-cell reference) and its use in the DLCO models.
Reference: GSE136831 (Adams/Sauler, Kaminski lab; same lab as GSE47460), Control + COPD cells only, 18 cell groups,
pseudobulk = mean over cells of counts/nUMI*1e4 (CP10K). The 2.0 GB count matrix was streamed and aggregated in the
user's browser pane (GEO is not reachable from the analysis container); the aggregated marker matrix is the input here.
Deconvolution: marker-based estimator of Hunt et al. (dtangle, Bioinformatics 2019) for log-scale arrays:
for each cell group k with specific markers M_k (top 20 by log2 specificity > 1 among genes on the array),
a_k = mean_{g in M_k}(y_g - log2(ref_gk + 1)), proportion p_k = 2^(a_k/gamma) / sum_j 2^(a_j/gamma), gamma = 1
(sensitivity gamma = 0.7). Probe-specific offsets make absolute proportions approximate; within-platform variation
across samples is what enters the models.
Outputs: estimated proportions per sample, validation against QC outliers and GSVA scores, and H5 re-run in R (10b).
Usage: python 10_deconvolution.py <prepared dir> <reference json> <cohort tsv (main 186)> <qc outlier list> <out dir>"""
import sys, os, json, gzip, numpy as np, pandas as pd
from scipy.optimize import nnls
prep, refj, cohort, qcf, out = sys.argv[1:6]; os.makedirs(out, exist_ok=True)
ref = json.load(open(refj)); groups = ref["groups"]
R = pd.DataFrame(ref["mean_cp10k"], index=ref["genes"], columns=groups)
R = R[~R.index.duplicated()]
coh = pd.read_csv(cohort, sep="\t"); qc = set(open(qcf).read().split())
rows = []; props = []
for p in ["GPL14550", "GPL6480"]:
    e = pd.read_csv(os.path.join(prep, f"{p}_expression.tsv.gz"), sep="\t", index_col=0)
    a = pd.read_csv(os.path.join(prep, f"{p}_gene_annotation.tsv.gz"), sep="\t", dtype=str, quoting=3)
    sym = a.set_index("ID")["GENE_SYMBOL"].reindex(e.index)
    s = coh[coh.platform == p]; E = e[s.sample_id]
    E = E.groupby(sym.values).mean()                       # one row per symbol (mean of probes)
    # marker selection: per group, genes present on the array with highest specificity (log2 ratio to max of others)
    common = R.index.intersection(E.index); Rc = R.loc[common]
    lr = np.log2(Rc + 1); mk = []; mks = {}
    for g in groups:
        spec = lr[g] - lr.drop(columns=g).max(axis=1)
        cand = spec[(spec > 1) & (Rc[g] >= 1)].sort_values(ascending=False).head(20).index.tolist(); mk += cand; mks[g] = cand
        rows.append(dict(platform=p, group=g, n_markers=len(cand)))
    mk = list(dict.fromkeys(mk))
    for gam in [1.0, 0.7]:
        A = np.vstack([(E.loc[mks[g]].values - np.log2(Rc.loc[mks[g], g].values[:, None] + 1)).mean(axis=0) for g in groups]).T
        W = 2 ** (A / gam); P = W / W.sum(axis=1, keepdims=True)
        if gam == 1.0:
            d = pd.DataFrame(P, columns=groups); d.insert(0, "sample_id", s.sample_id.values); d.insert(1, "platform", p)
            d["qc_outlier"] = d.sample_id.isin(qc)
            for j, g in enumerate(groups): d["logidx_" + g] = A[:, j]
        else:
            for j, g in enumerate(groups): d["g07_" + g] = P[:, j]
    props.append(d)
    pd.Series(mk).to_csv(os.path.join(out, f"markers_{p}.txt"), index=False, header=False)
D = pd.concat(props); D.to_csv(os.path.join(out, "deconv_proportions_186.csv"), index=False)
pd.DataFrame(rows).to_csv(os.path.join(out, "deconv_marker_counts.csv"), index=False)
D["alv_capillary_epi"] = D[["AT1", "AT2", "CapA", "CapB"]].sum(axis=1)
print(pd.DataFrame(rows).pivot(index="group", columns="platform", values="n_markers"))
print("\nMedian proportion (%) by platform, retained vs QC-excluded:")
print((D.groupby(["platform", "qc_outlier"])[groups + ["alv_capillary_epi"]].median().T * 100).round(1))
D["alv07"] = D[["g07_AT1", "g07_AT2", "g07_CapA", "g07_CapB"]].sum(axis=1)
print("\nSpearman gamma 1 vs 0.7, alveolar-capillary sum:", D.groupby("platform").apply(lambda x: x.alv_capillary_epi.corr(x.alv07, method="spearman")).round(3).to_dict())
