# GSE136831 reference aggregation (done in the browser, 2026-10-29)

The analysis container cannot reach GEO, so the 2.0 GB count matrix
`GSE136831_RawCounts_Sparse.mtx.gz` (45,947 genes x 312,928 cells, 692,789,348 non-zeros, all parsed)
was streamed from ftp.ncbi.nlm.nih.gov in 8 MB range requests and aggregated in the Claude desktop browser pane.

- Cells kept: Disease_Identity Control or COPD (IPF excluded), Manuscript_Identity mapped to 17 groups
  (Multiplet excluded; Macrophage_Alveolar + Macrophage merged).
  AT1=ATI; AT2=ATII; CapA=VE_Capillary_A (aerocyte; HPGD, EDNRB, SOSTDC1); CapB=VE_Capillary_B (general capillary; FCN3);
  VE_other=VE_Arterial/Venous/Peribronchial; Fibroblast=Fibroblast+Myofibroblast; SMC_Pericyte; Mesothelial; Lymphatic;
  Airway_epi=Basal/Club/Goblet/Ionocyte/PNEC/Aberrant_Basaloid; Ciliated; Monocyte=c/ncMonocyte;
  DC=cDC1/cDC2/pDC/DC_Mature/DC_Langerhans; Mast; T_NK=T/T_Cytotoxic/T_Regulatory/NK/ILC_A/ILC_B; B_Plasma.
- Pseudobulk per group = mean over cells of (count / nUMI x 1e4) (CP10K), nUMI from the metadata table.
- Exported genes: per group, up to 100 genes with log2(mean+1) - log2(max other+1) > 0.5 and mean >= 0.5 CP10K
  (union 1,404 genes) -> gse136831_reference_markers.json.
