#!/usr/bin/env bash
# Run the full GSE47460 COPD DLCO pipeline. Execute from the data folder (GEO_GSE47460_COPD/):
#   bash code/run_all.sh            (skip raw-file step 5 with: SKIP_RAW=1 bash code/run_all.sh)
set -euo pipefail
PY=${PY:-python3}
C=code
$PY VERIFY_FILES.py
(cd analysis/genesets && $PY ../../$C/python/make_gmt.py)
$PY $C/python/01_cohort_table1.py prepared/all_sample_metadata.tsv analysis/01_cohort
Rscript $C/R/02_h1_limma.R prepared analysis/01_cohort analysis/02_h1
Rscript $C/R/03_replication_diagnostics.R prepared analysis/01_cohort analysis/03_diagnostics
$PY $C/python/03b_array_qc_metrics.py prepared analysis/01_cohort analysis/03_diagnostics
Rscript $C/R/04_main_h1_to_h5.R prepared analysis/01_cohort analysis/genesets/travaglini_lung_c8.gmt analysis/04_main
if [ "${SKIP_RAW:-0}" != "1" ]; then
  mkdir -p expansion/RAW
  [ -n "$(ls -A expansion/RAW 2>/dev/null)" ] || tar -xf expansion/GSE47460_RAW.tar -C expansion/RAW
  $PY $C/python/05_raw_qc_batch.py expansion/RAW analysis/01_cohort/cohort_main_186.tsv analysis/03_diagnostics analysis/05_raw_qc
  Rscript $C/R/05b_batch_tissue.R prepared analysis/01_cohort analysis/genesets/travaglini_lung_c8.gmt analysis/05_raw_qc/raw_fe_headers_cohort186.csv analysis/05_raw_qc
fi
$PY $C/python/06_table1_qc177.py analysis/01_cohort analysis/03_diagnostics analysis/06_manuscript
$PY $C/python/07_manuscript_figures.py analysis analysis/03_diagnostics analysis/06_manuscript/figures
echo "done"
