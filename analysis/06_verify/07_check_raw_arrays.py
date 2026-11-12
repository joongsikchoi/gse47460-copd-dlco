"""Step 7: open and check every staged GSE47460 raw Agilent FE file.
Per file: gzip integrity, FEPARAMS/STATS/FEATURES blocks, platform grid, feature count, filename vs GEO title,
and identity: correlation of log2(gProcessedSignal) (non-control probes, mean per ProbeName) with every column of the
official series matrix of its platform -> best-matching column must be its own GSM.
Usage: python 07_check_raw_arrays.py <package dir> <device listing json> <out dir>"""
import sys, os, gzip, io, json, numpy as np, pandas as pd
pk, listing, out = sys.argv[1:4]
R = os.path.join(pk, "expansion/RAW")
expected = {x["name"]: x["size"] for x in json.load(open(listing))["entries"]}

def matrix(p):
    lines = gzip.open(os.path.join(pk, f"raw/GSE47460-{p}_series_matrix.txt.gz"), "rt").read().split("\n")
    b = lines.index("!series_matrix_table_begin"); e = lines.index("!series_matrix_table_end")
    title = [l for l in lines if l.startswith("!Sample_title")][0].split("\t")[1:]
    df = pd.read_csv(io.StringIO("\n".join(lines[b + 1:e])), sep="\t", index_col=0)
    df.columns = [c.strip('"') for c in df.columns]; df.index = [str(i).strip('"') for i in df.index]
    return df, dict(zip(df.columns, [t.strip('"') for t in title]))
M = {p: matrix(p) for p in ["GPL14550", "GPL6480"]}
grid = {"028004_D_F_20100430": "GPL14550", "014850_D_20070820": "GPL6480"}

rows = []
for fn in sorted(os.listdir(R)):
    gsm = fn.split("_")[0]; r = {"file": fn, "sample_id": gsm, "bytes": os.path.getsize(os.path.join(R, fn))}
    r["size_matches_device"] = r["bytes"] == expected.get(fn)
    try:
        text = gzip.open(os.path.join(R, fn), "rt", errors="replace").read(); r["gzip_ok"] = True
    except Exception as ex:
        r["gzip_ok"] = False; r["error"] = str(ex)[:80]; rows.append(r); continue
    blocks, hdr, feat = {}, None, []
    for line in text.split("\n"):
        p = line.split("\t")
        if p[0] in ("FEPARAMS", "STATS", "FEATURES"): hdr = (p[0], p[1:])
        elif p[0] == "DATA" and hdr:
            if hdr[0] == "FEATURES": feat.append(p[1:])
            elif hdr[0] not in blocks: blocks[hdr[0]] = dict(zip(hdr[1], p[1:]))
    fp = blocks.get("FEPARAMS", {})
    r.update(has_FEPARAMS="FEPARAMS" in blocks, has_STATS="STATS" in blocks, n_features=len(feat),
             grid=fp.get("Grid_Name"), scan_date=fp.get("Scan_Date", "")[:10], fe_version=fp.get("FeatureExtractor_Version"),
             barcode=fp.get("FeatureExtractor_Barcode"))
    plat = grid.get(r["grid"]); r["platform_from_grid"] = plat
    mx, titles = M[plat] if plat else (None, {})
    r["in_matrix_of_that_platform"] = gsm in (mx.columns if mx is not None else [])
    t = titles.get(gsm, ""); r["geo_title"] = t
    r["filename_matches_title"] = fn.replace(".txt.gz", "") == f"{gsm}_{t}"
    if feat and r["in_matrix_of_that_platform"]:
        f = pd.DataFrame(feat, columns=hdr[1])
        f = f[f.ControlType == "0"]; v = pd.to_numeric(f.gProcessedSignal, errors="coerce")
        s = np.log2(v.clip(lower=1)).groupby(f.ProbeName.values).mean()
        common = mx.index.intersection(s.index); r["probes_matched"] = len(common)
        c = mx.loc[common].corrwith(s.loc[common])
        r["cor_own"] = round(c[gsm], 4); best = c.drop(gsm).idxmax(); r["best_other"] = best; r["cor_best_other"] = round(c[best], 4)
        r["identity_ok"] = c.idxmax() == gsm
    rows.append(r)
d = pd.DataFrame(rows); d.to_csv(os.path.join(out, "raw_array_check.csv"), index=False)
print(len(d), "files checked")
for k in ["size_matches_device", "gzip_ok", "has_FEPARAMS", "has_STATS", "in_matrix_of_that_platform", "filename_matches_title", "identity_ok"]:
    print(f"{k:28s} {int(d[k].sum())}/{len(d)}")
print(d.groupby("platform_from_grid").agg(n=("file", "size"), feats_min=("n_features", "min"), feats_max=("n_features", "max"),
      cor_own_min=("cor_own", "min"), cor_own_median=("cor_own", "median"), gap_min=("cor_best_other", "max")))
print("FE versions:", d.fe_version.value_counts().to_dict())
bad = d[(~d.filename_matches_title) | (~d.identity_ok.fillna(False)) | (~d.gzip_ok)]
print(bad[["file", "geo_title", "cor_own", "best_other", "cor_best_other"]].to_string())
d["margin"] = d.cor_own - d.cor_best_other; print("smallest own-vs-other margins:"); print(d.nsmallest(8, "margin")[["file", "cor_own", "best_other", "cor_best_other", "margin"]].to_string())
