# Step 3: Why does the FEV1/emphysema-adjusted DLCO signal not replicate in GPL6480?
# (1) internal split-half replication within GPL14550 (test half n=61 = replication size)
# (2) power: expected direction concordance in GPL6480 given its SE, under full and shrunken effects
# (3) design-matched rerun (GOLD 1-4 only, + upper-lobe covariate, sex-discordant excluded)
# (4) hidden-factor adjustment (SVA) in each platform
# Usage: Rscript 03_replication_diagnostics.R <prepared dir> <01_cohort dir> <out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; out <- args[3]; dir.create(out, FALSE, TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
suppressPackageStartupMessages(library(sva)); set.seed(20261020)

sD <- main[main$platform == "GPL14550", ]; sR <- main[main$platform == "GPL6480", ]
ED <- P$GPL14550$e[genes, sD$sample_id]; ER <- P$GPL6480$e[genes, sR$sample_id]
D <- fitD(ED, sD); R <- fitD(ER, sR)
agree <- function(tA, pA, tB, k = 500) { top <- order(pA)[1:k]
  c(conc_top500 = mean(sign(tA[top]) == sign(tB[top])), tcor_top500 = cor(tA[top], tB[top]), tcor_all = cor(tA, tB)) }
cross <- agree(D$t, D$p, R$t)

# ---- (1) split-half within discovery ----
B <- 200
sh <- t(replicate(B, {
  idx <- sample(nrow(sD), nrow(sR)); a <- sD[-idx, ]; b <- sD[idx, ]
  fa <- fitD(ED[, a$sample_id], a); fb <- fitD(ED[, b$sample_id], b)
  c(nFDR05_train = sum(p.adjust(fa$p, "BH") < 0.05), agree(fa$t, fa$p, fb$t))
}))
sh_unadj <- t(replicate(100, {
  idx <- sample(nrow(sD), nrow(sR)); a <- sD[-idx, ]; b <- sD[idx, ]; f0 <- ~ DLCO + Age + Male + Current
  fa <- fitD(ED[, a$sample_id], a, f0); fb <- fitD(ED[, b$sample_id], b, f0); agree(fa$t, fa$p, fb$t)
}))
split_tab <- rbind(
  data.frame(comparison = "Cross-platform (GPL14550 -> GPL6480), adjusted", t(cross), q05 = NA, q95 = NA),
  data.frame(comparison = "Split-half within GPL14550 (64 -> 61), adjusted, median of 200", t(apply(sh[, -1], 2, median)),
             q05 = paste(round(apply(sh[, -1], 2, quantile, .05), 3), collapse = "/"), q95 = paste(round(apply(sh[, -1], 2, quantile, .95), 3), collapse = "/")),
  data.frame(comparison = "Split-half within GPL14550, unadjusted for FEV1/emph, median of 100", t(apply(sh_unadj, 2, median)), q05 = NA, q95 = NA))
split_tab$median_nFDR05_train_half <- c(NA, median(sh[, 1]), NA)
split_tab$frac_splits_conc_below_crossplatform <- c(NA, mean(sh[, "conc_top500"] <= cross["conc_top500"]), NA)
write.csv(split_tab, file.path(out, "diag1_split_half.csv"), row.names = FALSE); print(split_tab)

# ---- (2) power / expected concordance ----
top <- order(D$p)[1:500]; fdr <- which(p.adjust(D$p, "BH") < 0.05)
pw <- do.call(rbind, lapply(c(1, 0.5, 0.25), function(shrink) {
  zt <- abs(D$b * shrink) / R$se
  data.frame(effect_assumed = paste0(shrink * 100, "% of discovery estimate"),
    exp_conc_FDR05 = mean(pnorm(zt[fdr])), exp_conc_top500 = mean(pnorm(zt[top])),
    exp_power_p05_FDR05 = mean(pnorm(zt[fdr] - qnorm(.975)) + pnorm(-zt[fdr] - qnorm(.975))))
}))
se_ratio <- median(R$se / D$se)
pw$median_SE_ratio_rep_vs_disc <- se_ratio
pw$observed_conc_FDR05 <- mean(sign(D$t[fdr]) == sign(R$t[fdr])); pw$observed_conc_top500 <- cross["conc_top500"]
write.csv(pw, file.path(out, "diag2_power_expected_concordance.csv"), row.names = FALSE); print(pw)

# ---- (3) design-matched ----
FORM_L <- ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current + UpperLobe
dm <- function(filt, form) { a <- sD[filt(sD), ]; b <- sR[filt(sR), ]
  fa <- fitD(ED[, a$sample_id], a, form); fb <- fitD(ER[, b$sample_id], b, form)
  c(n_disc = nrow(a), n_rep = nrow(b), nFDR05_disc = sum(p.adjust(fa$p, "BH") < .05), nFDR05_rep = sum(p.adjust(fb$p, "BH") < .05), agree(fa$t, fa$p, fb$t)) }
okG <- function(s) is.na(s$GOLD) == FALSE & s$GOLD != "At Risk"
design_tab <- rbind(
  main_model = dm(function(s) rep(TRUE, nrow(s)), FORM),
  plus_upper_lobe = dm(function(s) rep(TRUE, nrow(s)), FORM_L),
  GOLD1to4_only = dm(okG, FORM),
  excl_sex_discordant = dm(function(s) !s$sex_discordant, FORM),
  all_three = dm(function(s) okG(s) & !s$sex_discordant, FORM_L),
  upper_lobe_samples_only = dm(function(s) s$UpperLobe == 1, FORM))
write.csv(design_tab, file.path(out, "diag3_design_matched.csv")); print(round(design_tab, 3))

# ---- (4) SVA ----
sv_fit <- function(E, s) { X <- model.matrix(FORM, s); X0 <- X[, colnames(X) != "DLCO"]
  n <- num.sv(E, X, method = "be", B = 20); sv <- if (n > 0) sva(E, X, X0, n.sv = n, B = 5)$sv else NULL
  s2 <- cbind(s, if (n > 0) setNames(as.data.frame(sv), paste0("SV", 1:n)) else NULL)
  form <- if (n > 0) update(FORM, as.formula(paste("~ . +", paste0("SV", 1:n, collapse = "+")))) else FORM
  r <- fitD(E, s2, form); r$nsv <- n; r }
invisible(capture.output(SD <- sv_fit(ED, sD), SR <- sv_fit(ER, sR)))
sva_tab <- data.frame(n_SV_disc = SD$nsv, n_SV_rep = SR$nsv, nFDR05_disc = sum(p.adjust(SD$p, "BH") < .05),
  nFDR05_rep = sum(p.adjust(SR$p, "BH") < .05), t(agree(SD$t, SD$p, SR$t)), tcor_disc_main_vs_sva = cor(D$t, SD$t))
write.csv(sva_tab, file.path(out, "diag4_sva.csv"), row.names = FALSE); print(sva_tab)

# ---- figure: split-half null distribution vs cross-platform ----
png(file.path(out, "fig4_split_half_vs_crossplatform.png"), width = 2400, height = 1100, res = 300)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1), las = 1, bty = "l")
for (m in c("conc_top500", "tcor_top500")) {
  h <- sh[, m]; hist(h, breaks = 25, col = "#2a78d6", border = "white", main = "", xlim = range(c(h, cross[m], if (m == "conc_top500") 0.4 else -0.2)),
    xlab = if (m == "conc_top500") "Direction concordance, top 500 genes" else "t correlation, top 500 genes", ylab = "Split-half replicates")
  abline(v = cross[m], col = "#eb6834", lwd = 2.5)
  title(main = if (m == "conc_top500") "Within GPL14550 (blue) vs GPL6480 (orange line)" else "", adj = 0, cex.main = 0.85, font.main = 1)
}
dev.off()

# ---- (6) split-half repeated on the array-QC cohort (fair comparison with QC primary analysis) ----
sDq <- sD[!sD$qc_outlier, ]; sRq <- sR[!sR$qc_outlier, ]
Dq <- fitD(ED[, sDq$sample_id], sDq); Rq <- fitD(ER[, sRq$sample_id], sRq); crossq <- agree(Dq$t, Dq$p, Rq$t)
shq <- t(replicate(200, { idx <- sample(nrow(sDq), nrow(sRq)); a <- sDq[-idx, ]; b <- sDq[idx, ]
  fa <- fitD(ED[, a$sample_id], a); fb <- fitD(ED[, b$sample_id], b); agree(fa$t, fa$p, fb$t) }))
qtab <- data.frame(comparison = c("Cross-platform, QC cohort", "Split-half within GPL14550 QC (59 -> 59), median of 200"),
  rbind(crossq, apply(shq, 2, median)), frac_splits_below_cross = c(NA, mean(shq[, "conc_top500"] <= crossq["conc_top500"])),
  conc_q05_q95 = c(NA, paste(round(quantile(shq[, "conc_top500"], c(.05, .95)), 3), collapse = "-")))
write.csv(qtab, file.path(out, "diag6_split_half_QC.csv"), row.names = FALSE); print(qtab)
