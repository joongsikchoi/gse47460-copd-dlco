# Step 14: cell number vs per-cell expression (exploratory, added after professor review 2026-11-08).
# For each alveolar-capillary group, take cell-enriched genes NOT used as deconvolution markers (no group's top-20 dtangle markers),
# score them (mean within-platform z), and test DLCO before and after adjusting for that group's estimated fraction.
# If the DLCO association of the non-marker score disappears after fraction adjustment -> fewer cells; if it persists -> lower expression per cell
# (or error in the fraction estimate). Main cohort QC177; all186 as sensitivity. Platforms fitted separately, IVW meta.
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; ddir <- args[3]; out <- args[4]; dir.create(out, FALSE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
suppressPackageStartupMessages(library(jsonlite))
ref <- fromJSON(file.path(ddir, "gse136831_reference_markers.json")); R <- ref$mean_cp10k; dimnames(R) <- list(ref$genes, ref$groups)
R <- R[!duplicated(rownames(R)), ]
dp <- read.csv(file.path(ddir, "deconv_proportions_186.csv")); fr <- c("AT1", "AT2", "CapA", "CapB")
dp[fr] <- 100 * dp[fr]; dp$AlvCap <- rowSums(dp[fr])
main <- merge(main, dp[, c("sample_id", setNames(fr, NULL), "AlvCap")], by = "sample_id", sort = FALSE)
markers <- unique(c(readLines(file.path(ddir, "markers_GPL14550.txt")), readLines(file.path(ddir, "markers_GPL6480.txt"))))
# specificity on the analysis gene universe (one probe per symbol)
com <- intersect(gene_sym, rownames(R)); lr <- log2(R[com, ] + 1)
sets <- list(); mk <- list()
for (g in fr) { spec <- lr[, g] - apply(lr[, colnames(lr) != g, drop = FALSE], 1, max)
  ok <- spec > 0.5 & R[com, g] >= 1 & !(com %in% markers)
  sets[[g]] <- com[ok][order(-spec[ok])][1:min(30, sum(ok))]                # up to 30 most enriched non-marker genes
  mk[[g]] <- intersect(markers, com[spec > 1 & R[com, g] >= 1]) }         # that group's own markers in the universe (positive control)
gs <- data.frame(group = rep(names(sets), lengths(sets)), gene = unlist(sets),
  spec_log2 = unlist(lapply(names(sets), function(g) (lr[sets[[g]], g] - apply(lr[sets[[g]], colnames(lr) != g, drop = FALSE], 1, max)))))
write.csv(gs, file.path(out, "nonmarker_gene_sets.csv"), row.names = FALSE); print(table(gs$group)); print(lengths(mk))
score <- function(p, s, sy) { id <- genes[match(sy, gene_sym)]; Z <- t(scale(t(P[[p]]$e[id, s$sample_id]))); colMeans(Z, na.rm = TRUE) }
models <- list(M0_base = "", M1_own_fraction = "+ FRAC", M2_AlvCap = "+ AlvCap", M3_four_fractions = "+ AT1 + AT2 + CapA + CapB")
meta <- function(b, se) { w <- 1 / se^2; B <- sum(w * b) / sum(w); S <- sqrt(1 / sum(w)); c(beta = B, se = S, p = 2 * pnorm(-abs(B / S))) }
res <- data.frame()
for (coh in c("QC177", "all186")) for (g in fr) for (typ in c("nonmarker", "marker")) for (m in names(models)) {
  sy <- if (typ == "nonmarker") sets[[g]] else mk[[g]]; if (!length(sy)) next
  st <- lapply(c("GPL14550", "GPL6480"), function(p) {
    s <- main[main$platform == p & (coh == "all186" | !main$qc_outlier), ]; s$y <- score(p, s, sy); s$FRAC <- s[[g]]
    f <- lm(as.formula(paste("y ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current", models[[m]])), s)
    c(b = 10 * coef(f)[["DLCO"]], se = 10 * sqrt(vcov(f)["DLCO", "DLCO"]), r_frac = cor(s$y, s$FRAC)) })
  mm <- meta(sapply(st, `[`, "b"), sapply(st, `[`, "se"))
  res <- rbind(res, data.frame(cohort = coh, group = g, genes = typ, n_genes = length(sy), model = m,
    beta_per10_disc = st[[1]][["b"]], beta_per10_rep = st[[2]][["b"]], meta_beta_per10 = mm[["beta"]], meta_se = mm[["se"]], meta_p = mm[["p"]],
    r_score_vs_fraction_disc = st[[1]][["r_frac"]], r_score_vs_fraction_rep = st[[2]][["r_frac"]])) }
res$pct_of_M0 <- with(res, 100 * meta_beta_per10 / ave(ifelse(model == "M0_base", meta_beta_per10, NA), cohort, group, genes, FUN = function(v) max(v, na.rm = TRUE)))
write.csv(res, file.path(out, "celltype_score_models.csv"), row.names = FALSE)
print(res[res$cohort == "QC177", c("group", "genes", "n_genes", "model", "meta_beta_per10", "meta_p", "pct_of_M0", "r_score_vs_fraction_disc")], digits = 3)
# gene level (non-marker genes, QC177): DLCO beta before vs after own-fraction adjustment
gl <- data.frame()
for (g in fr) { id <- genes[match(sets[[g]], gene_sym)]
  fx <- lapply(c("base", "frac"), function(v) lapply(c("GPL14550", "GPL6480"), function(p) {
    s <- main[main$platform == p & !main$qc_outlier, ]; s$FRAC <- s[[g]]; fitD(P[[p]]$e[id, s$sample_id], s, if (v == "base") FORM else update(FORM, ~ . + FRAC)) }))
  for (v in 1:2) { b <- sapply(fx[[v]], `[[`, "b"); se <- sapply(fx[[v]], `[[`, "se"); w <- 1 / se^2
    B <- rowSums(w * b) / rowSums(w); S <- sqrt(1 / rowSums(w)); fx[[v]] <- list(B = 10 * B, p = 2 * pnorm(-abs(B / S))) }
  gl <- rbind(gl, data.frame(group = g, gene = sets[[g]], beta_base = fx[[1]]$B, p_base = fx[[1]]$p, beta_fracadj = fx[[2]]$B, p_fracadj = fx[[2]]$p)) }
write.csv(gl, file.path(out, "nonmarker_gene_level.csv"), row.names = FALSE)
print(do.call(rbind, lapply(split(gl, gl$group), function(d) data.frame(group = d$group[1], n = nrow(d),
  pos_base = mean(d$beta_base > 0), nom_base = sum(d$p_base < .05), pos_adj = mean(d$beta_fracadj > 0), nom_adj = sum(d$p_fracadj < .05),
  median_ratio = median(d$beta_fracadj / d$beta_base)))), digits = 3)
