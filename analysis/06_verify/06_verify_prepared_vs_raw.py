"""Step 6: verify every prepared/ table against the official GEO raw files (series matrix, platform annotation).
Usage: python 06_verify_prepared_vs_raw.py <package dir> <out dir>"""
import sys, gzip, os, io, re, numpy as np, pandas as pd
pk, out = sys.argv[1], sys.argv[2]; R = []
def rep(check, ok, detail=""): R.append((check, "OK" if ok else "DIFF", detail)); print(("OK  " if ok else "DIFF"), check, detail)

def read_matrix(p):
    hdr, tab, intab = {}, [], False
    for line in gzip.open(p, "rt", encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if line.startswith("!series_matrix_table_begin"): intab = True; continue
        if line.startswith("!series_matrix_table_end"): break
        if intab: tab.append(line); continue
        if line.startswith("!"):
            k, *v = line.split("\t"); hdr.setdefault(k, []).append([x.strip('"') for x in v])
    df = pd.read_csv(io.StringIO("\n".join(tab)), sep="\t", index_col=0)
    df.columns = [c.strip('"') for c in df.columns]; df.index = [str(i).strip('"') for i in df.index]
    return hdr, df

allm = pd.read_csv(os.path.join(pk, "prepared/all_sample_metadata.tsv"), sep="\t", dtype=str, keep_default_na=False)
for p in ["GPL14550", "GPL6480"]:
    hdr, mx = read_matrix(os.path.join(pk, f"raw/GSE47460-{p}_series_matrix.txt.gz"))
    gsm = hdr["!Sample_geo_accession"][0]
    rep(f"{p} matrix: samples", True, f"{len(gsm)} samples, {mx.shape[0]} rows x {mx.shape[1]} cols, platform={set(hdr['!Sample_platform_id'][0])}")
    # expression
    ex = pd.read_csv(os.path.join(pk, f"prepared/{p}_expression.tsv.gz"), sep="\t", index_col=0)
    rep(f"{p} expression: same sample IDs & order", list(ex.columns) == gsm and list(mx.columns) == gsm)
    rep(f"{p} expression: same probe IDs & order", list(ex.index.astype(str)) == list(mx.index))
    mxv = mx.apply(pd.to_numeric, errors="coerce"); d = (ex.values - mxv.values)
    rep(f"{p} expression: values identical", np.nanmax(np.abs(d)) == 0 and np.array_equal(np.isnan(ex.values), np.isnan(mxv.values)),
        f"max |diff|={np.nanmax(np.abs(d)):.3g}, NA prepared={int(np.isnan(ex.values).sum())}, NA matrix={int(np.isnan(mxv.values).sum())}")
    # metadata: rebuild from matrix characteristics
    m = pd.DataFrame({"sample_id": gsm, "title": hdr["!Sample_title"][0], "source_name_ch1": hdr["!Sample_source_name_ch1"][0]})
    # characteristics sit in different row positions per sample -> parse each sample's key: value pairs
    per = [dict(x.split(": ", 1) for x in col if ": " in x) for col in zip(*hdr["!Sample_characteristics_ch1"])]
    for k in sorted({k for d in per for k in d}): m[k] = [d.get(k, "") for d in per]
    pm = pd.read_csv(os.path.join(pk, f"prepared/{p}_sample_metadata.tsv"), sep="\t", dtype=str, keep_default_na=False)
    rep(f"{p} sample metadata: same samples", list(pm.sample_id) == gsm)
    cols = [c for c in pm.columns if c not in ("sample_id", "platform")]
    miss = [c for c in cols if c not in m.columns]; extra = [c for c in m.columns if c not in pm.columns and c != "sample_id"]
    rep(f"{p} sample metadata: same columns as matrix characteristics", not miss and not extra, f"missing={miss} extra={extra}")
    diffs = {c: int((pm[c].str.strip() != m[c].str.strip()).sum()) for c in cols if c in m.columns}
    rep(f"{p} sample metadata: every value identical to matrix", sum(diffs.values()) == 0, str({k: v for k, v in diffs.items() if v}))
    # per-platform vs all_sample_metadata
    sub = allm[allm.platform == p].set_index("sample_id").loc[pm.sample_id]
    dd = {c: int((sub[c].values != pm[c].values).sum()) for c in pm.columns if c in sub.columns and c != "sample_id"}
    rep(f"{p}: per-platform metadata == all_sample_metadata rows", sum(dd.values()) == 0, str({k: v for k, v in dd.items() if v}))
    # annotation
    raw = gzip.open(os.path.join(pk, f"raw/{p}_annotation.txt.gz"), "rt", encoding="utf-8", errors="replace").read().split("\n")
    b = raw.index("!platform_table_begin"); e = raw.index("!platform_table_end")
    ra = pd.read_csv(io.StringIO("\n".join(raw[b + 1:e])), sep="\t", dtype=str, keep_default_na=False)
    pa = pd.read_csv(os.path.join(pk, f"prepared/{p}_gene_annotation.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
    rep(f"{p} annotation: rows", len(ra) == len(pa), f"raw {len(ra)} vs prepared {len(pa)}, columns equal={list(ra.columns) == list(pa.columns)}")
    same = ra.reset_index(drop=True).equals(pa.reset_index(drop=True))
    rep(f"{p} annotation: every cell identical to official platform table", same)
    rep(f"{p} annotation: all expression probes annotated", ex.index.astype(str).isin(pa.ID).all(), f"symbol missing for {int((pa.set_index('ID').loc[ex.index.astype(str), 'GENE_SYMBOL'] == '').sum())} of {len(ex)} probes")
    # series header vs source_series_metadata.txt
    if p == "GPL14550":
        ssm = open(os.path.join(pk, "source_series_metadata.txt"), encoding="utf-8").read()
        sids = re.findall(r"!Series_sample_id = (GSM\d+)", ssm)
        rep("source_series_metadata: sample IDs", len(sids) == 582 and set(sids) == set(allm.sample_id), f"{len(sids)} listed, {len(set(sids) & set(allm.sample_id))} match all_sample_metadata")
        for k in ["!Series_title", "!Series_summary", "!Series_overall_design"]:
            v = hdr[k][0][0]; rep(f"source_series_metadata: {k} matches matrix", v.strip() in ssm.replace("\n", " "))
rep("all_sample_metadata: 582 rows = 429 + 153", len(allm) == 582 and allm.platform.value_counts().to_dict() == {"GPL14550": 429, "GPL6480": 153})
pd.DataFrame(R, columns=["check", "result", "detail"]).to_csv(os.path.join(out, "verify_prepared_vs_raw.csv"), index=False)
