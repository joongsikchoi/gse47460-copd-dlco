# Step 4b: H2 platform estimates + inverse-variance meta for sensitivity cohorts (fills Table 3 "AT1 meta p" column).
# GSVA is re-scored within each cohort and platform (z-scored within platform), same model as 04_main.
# Usage: Rscript 04b_h2_meta_sensitivity.R <prepared dir> <01_cohort dir> <gmt> <out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; gmtf <- args[3]; out <- args[4]; dir.create(out, FALSE, TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
suppressPackageStartupMessages(library(GSVA))
cohorts <- list(
  primary_QC       = !main$qc_outlier,
  QC_GOLD1to4      = !main$qc_outlier & !(main$GOLD %in% "At Risk"),
  QC_sexconcordant = !main$qc_outlier & !main$sex_discordant,
  noQC_186         = rep(TRUE, nrow(main)))
gmt <- strsplit(readLines(gmtf), "\t"); gs <- setNames(lapply(gmt, `[`, -(1:2)), sub("\\|.*", "", sapply(gmt, `[`, 1)))
h2sets <- grep("^H2_", names(gs), value = TRUE)
res <- data.frame()
for (cn in names(cohorts)) {
  est <- list()
  for (p in names(P)) {
    s <- main[main$platform == p & cohorts[[cn]], ]; E <- P[[p]]$e[genes, s$sample_id]; rownames(E) <- gene_sym
    sc <- t(scale(t(gsva(gsvaParam(E, gs[h2sets], minSize = 10), verbose = FALSE)))); s <- cbind(s, t(sc))
    est[[p]] <- t(sapply(h2sets, function(g) summary(lm(reformulate(c("DLCO", "FEV1", "sqrtEmph", "Age", "Male", "Current"), g), s))$coef["DLCO", 1:2] * 10))
    attr(est[[p]], "n") <- nrow(s)
  }
  for (g in h2sets) {
    b <- c(est$GPL14550[g, 1], est$GPL6480[g, 1]); se <- c(est$GPL14550[g, 2], est$GPL6480[g, 2]); w <- 1 / se^2
    mb <- sum(w * b) / sum(w); ms <- sqrt(1 / sum(w)); Q <- sum(w * (b - mb)^2)
    res <- rbind(res, data.frame(cohort = cn, n_disc = attr(est$GPL14550, "n"), n_rep = attr(est$GPL6480, "n"), set = sub("H2_", "", g),
      beta_disc = b[1], se_disc = se[1], beta_rep = b[2], se_rep = se[2], meta_beta_per10 = mb, meta_se = ms,
      meta_p = 2 * pnorm(-abs(mb / ms)), heterogeneity_Q_p = pchisq(Q, 1, lower.tail = FALSE)))
  }
}
res$meta_p_holm <- ave(res$meta_p, res$cohort, FUN = function(x) p.adjust(x, "holm"))
write.csv(res, file.path(out, "H2_meta_sensitivity.csv"), row.names = FALSE); print(res, digits = 3)
