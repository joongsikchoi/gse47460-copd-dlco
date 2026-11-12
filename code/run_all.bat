@echo off
REM Run the full GSE47460 COPD DLCO pipeline on Windows. Execute from the data folder (GEO_GSE47460_COPD\):
REM   code\run_all.bat
REM Requires python and Rscript on PATH. Step 5 needs expansion\GSE47460_RAW.tar extracted to expansion\RAW\ first.
setlocal
chcp 65001 > nul
set C=code
python VERIFY_FILES.py || goto :err
pushd analysis\genesets && python ..\..\%C%\python\make_gmt.py && popd || goto :err
python %C%\python\01_cohort_table1.py prepared\all_sample_metadata.tsv analysis\01_cohort || goto :err
Rscript %C%\R\02_h1_limma.R prepared analysis\01_cohort analysis\02_h1 || goto :err
Rscript %C%\R\03_replication_diagnostics.R prepared analysis\01_cohort analysis\03_diagnostics || goto :err
python %C%\python\03b_array_qc_metrics.py prepared analysis\01_cohort analysis\03_diagnostics || goto :err
Rscript %C%\R\04_main_h1_to_h5.R prepared analysis\01_cohort analysis\genesets\travaglini_lung_c8.gmt analysis\04_main || goto :err
if exist expansion\RAW\*.txt.gz (
  python %C%\python\05_raw_qc_batch.py expansion\RAW analysis\01_cohort\cohort_main_186.tsv analysis\03_diagnostics analysis\05_raw_qc || goto :err
  Rscript %C%\R\05b_batch_tissue.R prepared analysis\01_cohort analysis\genesets\travaglini_lung_c8.gmt analysis\05_raw_qc\raw_fe_headers_cohort186.csv analysis\05_raw_qc || goto :err
) else (
  echo [skip] expansion\RAW not extracted - step 5 skipped
)
python %C%\python\06_table1_qc177.py analysis\01_cohort analysis\03_diagnostics analysis\06_manuscript || goto :err
python %C%\python\07_manuscript_figures.py analysis analysis\03_diagnostics analysis\06_manuscript\figures || goto :err
echo done
exit /b 0
:err
echo FAILED
exit /b 1
