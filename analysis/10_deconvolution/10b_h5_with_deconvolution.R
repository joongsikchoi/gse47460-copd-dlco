# Step 10b: H5 repeated with single-cell-reference deconvolution (GSE136831) instead of gene-set proxies,
# plus a direct test of estimated alveolar-capillary proportions against residual DLCO. QC cohort 118 / 59.
# Usage: Rscript 10b_h5_with_deconvolution.R <prepared dir> <01_cohort dir> <04_main out> <10 out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; o4 <- args[3]; out <- args[4]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
R4 <- readRDS(file.path(o4, "results.rds")); SC <- R4$SC
dc <- read.csv(file.path(out, "deconv_proportions_186.csv"))
dc$Lymphoid <- dc$T_NK + dc$B_Plasma
dc$Airway <- dc$Airway_epi + dc$Ciliated; dc$Stromal <- dc$Fibroblast + dc$SMC_Pericyte
dc$AlvCap <- dc$AT1 + dc$AT2 + dc$CapA + dc$CapB; dc$Cap <- dc$CapA + dc$CapB
comp <- c("Macrophage", "Lymphoid", "Airway", "Stromal")
h2sets <- c("H2_aCap", "H2_gCap", "H2_AT1", "H2_AT2"); clin <- c("FEV1", "sqrtEmph", "Age", "Male", "Current")
for (p in names(SC)) { s <- SC[[p]]; m <- dc[match(s$sample_id, dc$sample_id), ]; stopifnot(!anyNA(m$sample_id))
  for (v in c(comp, "AlvCap", "AT1", "AT2", "Cap", "CapA", "CapB")) s[[paste0("dc_", v)]] <- 100 * m[[v]]   # percent
  SC[[p]] <- s }
ivw <- function(b, se) { w <- 1 / se^2; m <- sum(w * b) / sum(w); s <- sqrt(1 / sum(w)); c(meta = m, se = s, p = 2 * pnorm(-abs(m / s)), Q_p = pchisq(sum(w * (b - m)^2), 1, lower.tail = FALSE)) }
cf <- function(s, y, covs) summary(lm(reformulate(c("DLCO", covs), y), s))$coef["DLCO", c(1, 2, 4)]
dcomp <- paste0("dc_", comp)

# (1) validation: deconvolution vs GSVA scores
val <- do.call(rbind, lapply(names(SC), function(p) { s <- SC[[p]]
  data.frame(platform = p, AT1_vs_GSVA_AT1 = cor(s$dc_AT1, s$H2_AT1, method = "s"), AT2_vs_GSVA_AT2 = cor(s$dc_AT2, s$H2_AT2, method = "s"),
             Cap_vs_GSVA_gCap = cor(s$dc_Cap, s$H2_gCap, method = "s"), Cap_vs_GSVA_aCap = cor(s$dc_Cap, s$H2_aCap, method = "s"),
             Mac_vs_GSVA_Mac = cor(s$dc_Macrophage, s$Macrophage, method = "s"), Airway_vs_GSVA_Cil = cor(s$dc_Airway, s$Ciliated, method = "s")) }))
write.csv(val, file.path(out, "deconv_validation_vs_gsva.csv"), row.names = FALSE); print(val, digits = 2)

# (2) H2 scores adjusted for deconvolved non-alveolar composition
h2 <- data.frame()
for (g in h2sets) { r <- sapply(names(SC), function(p) c(cf(SC[[p]], g, clin), cf(SC[[p]], g, c(clin, dcomp))))
  m0 <- ivw(r[1, ] * 10, r[2, ] * 10); m1 <- ivw(r[4, ] * 10, r[5, ] * 10)
  h2 <- rbind(h2, data.frame(set = sub("H2_", "", g), disc_beta = r[1, 1] * 10, disc_beta_adj = r[4, 1] * 10, disc_p_adj = r[6, 1],
    rep_beta = r[1, 2] * 10, rep_beta_adj = r[4, 2] * 10, rep_p_adj = r[6, 2], meta_beta = m0["meta"], meta_p = m0["p"],
    meta_beta_adj = m1["meta"], meta_p_adj = m1["p"], ratio_disc = r[4, 1] / r[1, 1], ratio_rep = r[4, 2] / r[1, 2])) }
write.csv(h2, file.path(out, "H5_deconv_H2.csv"), row.names = FALSE); print(h2, digits = 3)

# (3) gene level: H1 and top-100 discovery genes after deconvolution adjustment
D <- R4$h1; topD <- NULL
gl <- data.frame()
base <- list(); adj <- list()
for (p in names(SC)) { s <- SC[[p]]; E <- P[[p]]$e[genes, s$sample_id]
  base[[p]] <- fitD(E, s); adj[[p]] <- fitD(E, s, update(FORM, as.formula(paste("~ . +", paste(dcomp, collapse = "+"))))) }
top <- order(base$GPL14550$p)[1:100]
for (p in names(SC)) gl <- rbind(gl, data.frame(platform = p, FDR05_after = sum(p.adjust(adj[[p]]$p, "BH") < .05),
  tcor_vs_unadjusted = cor(base[[p]]$t, adj[[p]]$t), top100_same_direction = mean(sign(adj[[p]]$b[top]) == sign(base[[p]]$b[top])),
  top100_median_beta_ratio = median(adj[[p]]$b[top] / base[[p]]$b[top])))
topA <- order(adj$GPL14550$p)[1:100]; gl$conc_top100_adjusted_disc_to_rep <- mean(sign(adj$GPL14550$t[topA]) == sign(adj$GPL6480$t[topA]))
write.csv(gl, file.path(out, "H5_deconv_genes.csv"), row.names = FALSE); print(gl, digits = 3)

# (4) direct test: estimated proportions (%) ~ DLCO + clinical covariates
pr <- data.frame()
for (v in c("AlvCap", "AT1", "AT2", "Cap", "CapA", "CapB", comp)) { y <- paste0("dc_", v)
  r <- sapply(names(SC), function(p) cf(SC[[p]], y, clin)); m <- ivw(r[1, ] * 10, r[2, ] * 10)
  pr <- rbind(pr, data.frame(proportion = v, disc_median_pct = median(SC$GPL14550[[y]]), rep_median_pct = median(SC$GPL6480[[y]]),
    disc_beta_pct_per10 = r[1, 1] * 10, disc_p = r[3, 1], rep_beta_pct_per10 = r[1, 2] * 10, rep_p = r[3, 2], meta_beta = m["meta"], meta_p = m["p"], Q_p = m["Q_p"])) }
write.csv(pr, file.path(out, "deconv_proportions_vs_DLCO.csv"), row.names = FALSE); print(pr, digits = 3)
