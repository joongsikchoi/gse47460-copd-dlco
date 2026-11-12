# Step 13: post hoc diagnostics requested in the 2026-11-08 review (outlier handling, H5, multiplicity, effect size).
#  A. 186-patient cohort (no outlier exclusion) with the estimated alveolar-capillary fraction as covariate
#  B. limma arrayWeights in the 186-patient cohort
#  C. Cook's distance of each sample for the 375 discovery genes
#  D. alternative exclusion cutoffs (Q1 - k x IQR, k = 1.5, 2, 3)
#  E. influence of individual outlier arrays (leave-one-in / leave-subset-out)
#  F. prespecified H5 (non-alveolar GSVA scores) with inverse-variance meta-analysis
#  G. multiplicity (BH) across cell groups; effect size relative to the mean alveolar-capillary fraction
# Usage: Rscript 13_reviewer_diagnostics.R <prepared dir> <01_cohort dir> <gmt> <10_deconvolution dir> <out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; gmtf <- args[3]; ddir <- args[4]; out <- args[5]; dir.create(out, FALSE, TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
suppressPackageStartupMessages(library(GSVA)); set.seed(20261020)

dp <- read.csv(file.path(ddir, "deconv_proportions_186.csv"))
grp <- c("AT1","AT2","CapA","CapB","VE_other","Lymphatic","Fibroblast","SMC_Pericyte","Mesothelial","Airway_epi","Ciliated","Monocyte","DC","Mast","T_NK","B_Plasma","Macrophage")
dp[, grp] <- dp[, grp] * 100
dp$AlvCap <- dp$AT1 + dp$AT2 + dp$CapA + dp$CapB
dp <- dp[, c("sample_id", grp, "AlvCap")]; names(dp)[-1] <- paste0(names(dp)[-1], "_dt"); grp_dt <- paste0(grp, "_dt")
main <- merge(main, dp, by = "sample_id", sort = FALSE); main$AlvCap <- main$AlvCap_dt
get <- function(filt, p) main[main$platform == p & filt, ]
cnt <- function(f) sum(p.adjust(f$p, "BH") < .05)
conc <- function(D, R, k = 500) { s <- order(D$p)[1:k]; mean(sign(D$t[s]) == sign(R$t[s])) }
res <- list()

# ---------- A + B: 186 cohort, covariate / arrayWeights ----------
FA <- update(FORM, ~ . + AlvCap)
A <- data.frame()
for (p in names(P)) {
  s <- get(rep(TRUE, nrow(main)), p); E <- P[[p]]$e[genes, s$sample_id]
  base <- fitD(E, s); adj <- fitD(E, s, FA)
  X <- model.matrix(FORM, s); aw <- arrayWeights(E, X)
  fw <- eBayes(lmFit(E, X, weights = aw), trend = TRUE, robust = TRUE)
  A <- rbind(A, data.frame(platform = p, n = nrow(s), FDR05_base = cnt(base), FDR05_AlvCap_adj = cnt(adj),
    FDR05_arrayWeights = sum(p.adjust(fw$p.value[, "DLCO"], "BH") < .05),
    aw_outliers_median = median(aw[s$qc_outlier]), aw_others_median = median(aw[!s$qc_outlier]),
    aw_outliers_rank_max = max(rank(aw)[s$qc_outlier])))
  res[[p]] <- list(s = s, E = E, base = base, adj = adj, fw = fw, aw = setNames(aw, s$sample_id))
}
# cross-platform concordance of the 186 analyses
A$conc_top500_base <- conc(res$GPL14550$base, res$GPL6480$base)
A$conc_top500_AlvCap <- conc(res$GPL14550$adj, res$GPL6480$adj)
fwD <- res$GPL14550$fw; fwR <- res$GPL6480$fw
A$conc_top500_arrayWeights <- { s <- order(fwD$p.value[, "DLCO"])[1:500]; mean(sign(fwD$t[s, "DLCO"]) == sign(fwR$t[s, "DLCO"])) }
write.csv(A, file.path(out, "A_B_186_covariate_arrayweights.csv"), row.names = FALSE); print(A)
aw_tab <- data.frame(sample_id = names(res$GPL14550$aw), weight = res$GPL14550$aw)
aw_tab$outlier <- aw_tab$sample_id %in% main$sample_id[main$qc_outlier]
write.csv(aw_tab[order(aw_tab$weight), ], file.path(out, "B_arrayweights_GPL14550.csv"), row.names = FALSE)

# ---------- C: Cook's distance for the 375 discovery genes ----------
sD <- res$GPL14550$s; ED <- res$GPL14550$E; sel <- which(p.adjust(res$GPL14550$base$p, "BH") < .05)
X <- model.matrix(FORM, sD); n <- nrow(X)
CD <- sapply(sel, function(i) cooks.distance(lm(ED[i, ] ~ X - 1)))
rownames(CD) <- sD$sample_id
C <- data.frame(sample_id = sD$sample_id, outlier = sD$qc_outlier, DLCO = sD$DLCO,
  resid_DLCO = resid(lm(DLCO ~ FEV1 + sqrtEmph + Age + Male + Current, sD)), AlvCap = sD$AlvCap,
  median_cooksD = apply(CD, 1, median), pct_genes_cooksD_gt_4n = rowMeans(CD > 4 / n) * 100,
  pct_genes_cooksD_gt_1 = rowMeans(CD > 1) * 100)
C <- C[order(-C$median_cooksD), ]; write.csv(C, file.path(out, "C_cooksD_375genes.csv"), row.names = FALSE)
cat("n genes:", length(sel), " 4/n =", round(4 / n, 4), "\n"); print(head(C, 12), digits = 3)

# ---------- D: alternative cutoffs ----------
mcor <- function(p) { s <- main[main$platform == p, ]; E <- P[[p]]$e[genes, s$sample_id]; cc <- cor(E)
  setNames(apply(cc, 1, function(x) median(x[x < 1])), s$sample_id) }
MC <- c(mcor("GPL14550"), mcor("GPL6480"))
Dt <- data.frame()
for (k in c(1.5, 2, 3)) {
  ex <- unlist(lapply(c("GPL14550", "GPL6480"), function(p) { m <- MC[main$sample_id[main$platform == p]]; q <- quantile(m, c(.25, .75)); names(m)[m < q[1] - k * diff(q)] }))
  f <- !(main$sample_id %in% ex); sD2 <- get(f, "GPL14550"); sR2 <- get(f, "GPL6480")
  D2 <- fitD(P$GPL14550$e[genes, sD2$sample_id], sD2); R2 <- fitD(P$GPL6480$e[genes, sR2$sample_id], sR2)
  Dt <- rbind(Dt, data.frame(k_IQR = k, excl_disc = sum(main$platform[main$sample_id %in% ex] == "GPL14550"),
    excl_rep = sum(main$platform[main$sample_id %in% ex] == "GPL6480"), n_disc = nrow(sD2), n_rep = nrow(sR2),
    FDR05_disc = cnt(D2), FDR05_rep = cnt(R2), conc_top500 = conc(D2, R2)))
}
write.csv(Dt, file.path(out, "D_cutoffs.csv"), row.names = FALSE); print(Dt)

# ---------- E: influence of individual outlier arrays (discovery) ----------
outD <- sD$sample_id[sD$qc_outlier]; E1 <- data.frame()
for (o in outD) { keep <- !(sD$sample_id %in% setdiff(outD, o))          # 118 non-outliers + this one outlier
  s2 <- sD[keep, ]; E1 <- rbind(E1, data.frame(added_back = o, n = nrow(s2), FDR05 = cnt(fitD(ED[, s2$sample_id], s2)))) }
for (o in outD) { keep <- sD$sample_id != o                               # 125 minus this one outlier
  s2 <- sD[keep, ]; E1 <- rbind(E1, data.frame(added_back = paste0("all_but_", o), n = nrow(s2), FDR05 = cnt(fitD(ED[, s2$sample_id], s2)))) }
top2 <- head(C$sample_id[C$outlier], 2)
for (lab in c("drop_top2_only", "drop_other5_only")) {
  dr <- if (lab == "drop_top2_only") top2 else setdiff(outD, top2); s2 <- sD[!(sD$sample_id %in% dr), ]
  E1 <- rbind(E1, data.frame(added_back = paste(lab, paste(dr, collapse = ";")), n = nrow(s2), FDR05 = cnt(fitD(ED[, s2$sample_id], s2)))) }
write.csv(E1, file.path(out, "E_influence_outliers.csv"), row.names = FALSE); print(E1)

# ---------- F: prespecified H5 (GSVA proxies) with meta-analysis ----------
gmt <- strsplit(readLines(gmtf), "\t"); gs <- setNames(lapply(gmt, `[`, -(1:2)), sub("\\|.*", "", sapply(gmt, `[`, 1)))
h2sets <- grep("^H2_", names(gs), value = TRUE); comp4 <- c("Macrophage", "AlvFibroblast", "Neutrophil", "Ciliated")
Fr <- data.frame()
for (p in names(P)) { s <- get(!main$qc_outlier, p); Ex <- P[[p]]$e[genes, s$sample_id]; rownames(Ex) <- gene_sym
  sc <- t(scale(t(gsva(gsvaParam(Ex, gs, minSize = 10), verbose = FALSE)))); s <- cbind(s, t(sc))
  for (g in h2sets) for (adj in c("none", "GSVA_comp4", "dtangle_nonalv")) {
    cv <- switch(adj, none = NULL, GSVA_comp4 = comp4, dtangle_nonalv = c("Macrophage_dt", "Lymphoid_d", "Airway_d", "Stromal_d"))
    s$Lymphoid_d <- s$T_NK_dt + s$B_Plasma_dt; s$Airway_d <- s$Airway_epi_dt + s$Ciliated_dt; s$Stromal_d <- s$Fibroblast_dt + s$SMC_Pericyte_dt
    f <- summary(lm(reformulate(c("DLCO", "FEV1", "sqrtEmph", "Age", "Male", "Current", cv), g), s))$coef["DLCO", ]
    Fr <- rbind(Fr, data.frame(platform = p, set = sub("H2_", "", g), adjustment = adj, beta10 = f[1] * 10, se10 = f[2] * 10, p = f[4])) } }
meta <- do.call(rbind, lapply(split(Fr, list(Fr$set, Fr$adjustment)), function(d) { w <- 1 / d$se10^2; b <- sum(w * d$beta10) / sum(w); se <- sqrt(1 / sum(w))
  data.frame(set = d$set[1], adjustment = d$adjustment[1], disc_beta = d$beta10[d$platform == "GPL14550"], disc_p = d$p[d$platform == "GPL14550"],
    rep_beta = d$beta10[d$platform == "GPL6480"], rep_p = d$p[d$platform == "GPL6480"], meta_beta = b, meta_lo = b - 1.96 * se, meta_hi = b + 1.96 * se, meta_p = 2 * pnorm(-abs(b / se))) }))
meta <- meta[order(meta$set, meta$adjustment), ]; write.csv(meta, file.path(out, "F_H5_prespecified_vs_dtangle.csv"), row.names = FALSE); print(meta, digits = 3)

# ---------- G: multiplicity across cell groups + relative effect size ----------
G <- data.frame(); qc <- !main$qc_outlier
for (g in c(grp_dt, "AlvCap")) { r <- lapply(names(P), function(p) { s <- get(qc, p); summary(lm(reformulate(c("DLCO", "FEV1", "sqrtEmph", "Age", "Male", "Current"), g), s))$coef["DLCO", ] })
  b <- sapply(r, `[`, 1) * 10; se <- sapply(r, `[`, 2) * 10; w <- 1 / se^2; mb <- sum(w * b) / sum(w); ms <- sqrt(1 / sum(w))
  G <- rbind(G, data.frame(group = g, disc_beta = b[1], rep_beta = b[2], rep_p = sapply(r, `[`, 4)[2], meta_beta = mb, meta_p = 2 * pnorm(-abs(mb / ms)),
    mean_pct = mean(main[[g]][qc]), rel_change_per10_pct = 100 * mb / mean(main[[g]][qc]))) }
G$BH_17groups <- NA; i17 <- G$group %in% grp_dt; G$BH_17groups[i17] <- p.adjust(G$meta_p[i17], "BH")
write.csv(G, file.path(out, "G_cellgroups_multiplicity.csv"), row.names = FALSE); print(G, digits = 3)
cat("AlvCap mean % (QC 177):", round(mean(main$AlvCap[qc]), 2), " SD:", round(sd(main$AlvCap[qc]), 2), " DLCO SD:", round(sd(main$DLCO[qc]), 2), "\n")
