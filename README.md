# Alveolar-capillary cell composition and residual diffusing capacity in COPD lung tissue

Analysis code and results for the study of diffusing capacity (DLCO) that is **not** explained by
airflow limitation and CT emphysema in COPD lung tissue, using the public microarray dataset
**GSE47460** (Lung Genomics Research Consortium / Lung Tissue Research Consortium).

> **Status:** the manuscript is under preparation / review. This repository is released so that the
> analysis can be inspected and reproduced. The citation below will be updated on acceptance.

## What this study does

Among COPD patients with similar FEV1 and CT emphysema, some have a markedly lower DLCO. We asked
whether that residual variation is associated with (1) the expression of individual genes or
(2) the amount of gas-exchanging cells in the tissue, analysing the dataset's two microarray
platforms separately so that one serves as an internal replication of the other.

| | |
|---|---|
| Cohort | 220 COPD patients → 186 complete cases → **177** after excluding 9 tissue-composition outlier arrays |
| Platforms | GPL14550 (discovery, n = 118) and GPL6480 (replication, n = 59), analysed separately and combined by inverse-variance meta-analysis |
| Model | `Y ~ DLCO + FEV1 + sqrt(emphysema) + age + sex + current smoking`, fitted per platform; effects are per 10-point higher DLCO |
| Hypotheses | H1 gene-level association · H2 prespecified cell-type scores · H3 cross-platform concordance · H4 added predictive value · H5 robustness to cell composition |
| Headline result | No gene-level association replicated across platforms; the estimated alveolar-capillary cell fraction was 0.96 percentage points higher per 10-point higher DLCO (meta p = 0.0013) |

## Data

This repository contains **code and derived results only**. The primary data are public and are not
redistributed here:

| Dataset | Use | Where to get it |
|---|---|---|
| [GSE47460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460) | Lung tissue expression and clinical variables | NCBI GEO |
| [GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831) | Single-cell reference for deconvolution | NCBI GEO |

The clinical tables under `analysis/01_cohort/` are derived from the GEO sample metadata of
GSE47460, which is itself public and de-identified.

## Layout

```
code/          curated scripts and run_all.sh
  R/           limma models, hypothesis tests, diagnostics
  python/      cohort tables, array QC, figures
analysis/      one directory per analysis step, with the result files it produced
  01_cohort/   cohort selection, Table 1, clinical model
  02_h1/       first-pass gene-level analysis (before outlier exclusion)
  03_diagnostics/  why the first pass did not replicate; array QC rule
  04_main/     main analysis of H1-H5 in the 177-patient cohort
  05_raw_qc/   raw Agilent Feature Extraction files, scan-batch checks
  06_verify/   checksum and value-level verification of the prepared tables
  08_figures/  supplementary figures, split-half draws
  09_additional/   model variants, leave-one-out, camera
  10_deconvolution/  GSE136831 reference and dtangle cell fractions
  11_figure6/  cell-composition figure
  12_figures_english/  English-labelled versions of the figures
  13_reviewer/ post hoc diagnostics of the outlier rule, robust regression
  14_celltype_expression/  cell number versus expression per cell
  genesets/    prespecified MSigDB C8 gene sets (Travaglini et al.)
figures/       final manuscript figures (vector PDF)
```

## Reproducing the analysis

1. Download GSE47460 (series matrix and, for the quality checks, `GSE47460_RAW.tar`) and GSE136831
   from GEO, and place the normalised expression tables and annotation under a `prepared/` directory.
2. Install the dependencies:
   - **R** (≥ 4.3): `limma`, `GSVA`, `sva`, `MASS`, `org.Hs.eg.db`, `GO.db`
   - **Python** (≥ 3.10): `numpy`, `pandas`, `scipy`, `matplotlib`
3. Run the steps in order. `code/run_all.sh` runs the main pipeline (cohort table, first-pass gene
   analysis, replication diagnostics, main H1-H5, raw-file QC, Table 1 and the manuscript figures)
   and shows the arguments each script takes. The later steps, under `analysis/08_figures` through
   `analysis/14_celltype_expression`, are run individually and follow the same convention: inputs
   first, output directory last. Every script writes its outputs to that directory, which is what
   the `analysis/` subdirectories here contain.

The random seed is `20261020` throughout.

## Notes

- **Tissue-composition outliers.** The rule that excludes 9 arrays (median inter-array correlation
  below Q1 − 3×IQR) was introduced *after* the first-pass analysis, while diagnosing why the
  discovery signals did not appear on the other platform. It is reported as a post hoc diagnostic,
  and every main result is also reported without exclusion (`analysis/13_reviewer/`).
- **Figures.** `code/python/make_figures_hires.py` re-renders the figures at 600 dpi (TIFF) and as
  vector EPS/PDF. It expects `code/python/07_manuscript_figures_patched.py`, which differs from
  `07_manuscript_figures.py` only in deriving the outlier flag from the sample list rather than from
  a `qc_outlier` column that is absent from the saved QC tables.
- The EPS files that matplotlib produces render semi-transparent points as opaque, because the
  PostScript backend has no transparency; use the TIFF or PDF versions.

## Citation

> Choi JS, Kim K. Alveolar-capillary cell composition and diffusing capacity unexplained by airflow
> limitation and emphysema in COPD lung tissue. Manuscript under review.

## License

[MIT](LICENSE) © 2026 Joong Sik Choi, Kangjin Kim
