"""Step 5: parse Agilent Feature Extraction headers (FEPARAMS + STATS) of raw GSE47460 files for the COPD cohort.
Usage: python 05_raw_qc_batch.py <RAW dir> <cohort tsv> <qc_array csv dir> <out dir>"""
import sys, gzip, glob, os, re, pandas as pd
raw, coh, qcd, out = sys.argv[1:5]
rows = []
for f in sorted(glob.glob(os.path.join(raw, "*.txt.gz"))):
    gsm = os.path.basename(f).split("_")[0]; blocks = {}; hdr = None
    with gzip.open(f, "rt", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] in ("FEPARAMS", "STATS"): hdr = (p[0], p[1:])
            elif p[0] == "DATA" and hdr and hdr[0] not in blocks: blocks[hdr[0]] = dict(zip(hdr[1], p[1:])); 
            elif p[0] == "FEATURES": break
    r = {"sample_id": gsm, "file": os.path.basename(f)}
    fp = blocks.get("FEPARAMS", {}); r.update({k: fp.get(k) for k in ["Scan_Date", "Scan_ScannerName", "FeatureExtractor_Barcode", "FeatureExtractor_ScanFileName", "Grid_Name", "Protocol_Name", "FeatureExtractor_Version"]})
    r.update({"STATS_" + k: v for k, v in blocks.get("STATS", {}).items()})
    rows.append(r)
d = pd.DataFrame(rows)
d["scan_date"] = pd.to_datetime(d.Scan_Date.str.slice(0, 10), format="%m-%d-%Y", errors="coerce").dt.date
d["slide"] = d.FeatureExtractor_Barcode.str.split("_").str[0]
d["batch_folder"] = d.FeatureExtractor_ScanFileName.str.extract(r"\\\\?([^\\]+)\\\\?[^\\]*$")[0]
c = pd.read_csv(coh, sep="\t"); q = pd.concat([pd.read_csv(os.path.join(qcd, f"qc_array_{p}.csv")) for p in ["GPL14550", "GPL6480"]])
d = d.merge(c[["sample_id", "platform", "DLCO", "FEV1", "Emphysema", "GOLD"]], on="sample_id").merge(q[["sample_id", "medcor", "z_medcor"]], on="sample_id")
d.to_csv(os.path.join(out, "raw_fe_headers_cohort186.csv"), index=False)
print(d.shape); print(d.groupby("platform").agg(n=("sample_id", "size"), dates=("scan_date", "nunique"), slides=("slide", "nunique"), first=("scan_date", "min"), last=("scan_date", "max")))
print(d.Grid_Name.value_counts()); print(d.FeatureExtractor_Version.value_counts()); print(d.batch_folder.value_counts().head(20))
