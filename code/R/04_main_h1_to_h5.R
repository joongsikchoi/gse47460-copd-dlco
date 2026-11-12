# Step 4: Primary analysis H1-H5 after array QC.
# Primary cohort: main 186 minus array-QC outliers (Tukey far-out on median inter-array correlation, same rule on both platforms)
#   -> GPL14550 118 / GPL6480 59. QC rule was added after the first-pass analysis (step 2) and is disclosed as such.
# Gene sets: MSigDB C8 Travaglini 2020 (genesets/travaglini_lung_c8.gmt), fixed before this script was run.
# Usage: Rscript 04_main_h1_to_h5.R <prepared dir> <01_cohort dir> <gmt> <out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; gmtf <- args[3]; out <- args[4]; dir.create(out, FALSE, TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
suppressPackageStartupMessages(library(GSVA)); set.seed(20261020)

cohorts <- list(
  primary_QC        = !main$qc_outlier,
  noQC_186          = rep(TRUE, nrow(main)),
  QC_GOLD1to4       = !main$qc_outlier & !(main$GOLD %in% "At Risk"),
  QC_sexconcordant  = !main$qc_outlier & !main$sex_discordant)
get <- function(filt, p) main[main$platform == p & filt, ]

# ================= H1 + H3 =================
perm_rep <- function(Dres, sR, ER, sel, B = 1000) {           # Freedman-Lane permutation of DLCO in replication
  X <- model.matrix(FORM, sR); Z <- X[, colnames(X) != "DLCO"]
  h <- fitted(lm(X[, "DLCO"] ~ Z - 1)); e <- X[, "DLCO"] - h
  ds <- sign(Dres$t[sel]); E <- ER[sel, , drop = FALSE]
  st <- function(tv) c(mean(sign(tv) == ds), cor(Dres$t[sel], tv))
  obs <- st(fitD(E, sR)$t)
  nul <- t(replicate(B, { s2 <- sR; s2$DLCO <- h + sample(e); st(fitD(E, s2)$t) }))
  c(n = length(sel), concordance = obs[1], p_conc = (1 + sum(nul[, 1] >= obs[1])) / (B + 1),
    t_cor = obs[2], p_tcor = (1 + sum(nul[, 2] >= obs[2])) / (B + 1))
}
h1 <- list(); h1sum <- data.frame(); h3 <- data.frame()
for (cn in names(cohorts)) {
  sD <- get(cohorts[[cn]], "GPL14550"); sR <- get(cohorts[[cn]], "GPL6480")
  ED <- P$GPL14550$e[genes, sD$sample_id]; ER <- P$GPL6480$e[genes, sR$sample_id]
  D <- fitD(ED, sD); R <- fitD(ER, sR); qD <- p.adjust(D$p, "BH"); qR <- p.adjust(R$p, "BH")
  h1[[cn]] <- list(D = D, R = R, sD = sD, sR = sR)
  h1sum <- rbind(h1sum, data.frame(cohort = cn, n_disc = nrow(sD), n_rep = nrow(sR),
    disc_FDR05 = sum(qD < .05), disc_FDR10 = sum(qD < .10), disc_FDR25 = sum(qD < .25), disc_pi0 = round(propTrueNull(D$p), 3),
    rep_FDR05 = sum(qR < .05), rep_pi0 = round(propTrueNull(R$p), 3)))
  for (k in c(100, 500)) { sel <- order(D$p)[1:k]
    h3 <- rbind(h3, data.frame(cohort = cn, gene_set = paste0("top", k), t(perm_rep(D, sR, ER, sel, if (cn == "primary_QC") 1000 else 300)))) }
  if (sum(qD < .05) >= 10) h3 <- rbind(h3, data.frame(cohort = cn, gene_set = "FDR05", t(perm_rep(D, sR, ER, which(qD < .05), 300))))
}
write.csv(h1sum, file.path(out, "H1_summary.csv"), row.names = FALSE); write.csv(h3, file.path(out, "H3_replication.csv"), row.names = FALSE)
print(h1sum); print(h3, digits = 3)
D <- h1$primary_QC$D; R <- h1$primary_QC$R
tab <- data.frame(probe = genes, symbol = gene_sym, disc_beta_per10 = D$b * 10, disc_t = D$t, disc_p = D$p, disc_FDR = p.adjust(D$p, "BH"),
  rep_beta_per10 = R$b * 10, rep_t = R$t, rep_p = R$p, rep_FDR = p.adjust(R$p, "BH"),
  same_direction = sign(D$t) == sign(R$t), meta_z = (D$t * sqrt(118) + R$t * sqrt(59)) / sqrt(118 + 59))
tab$meta_p <- 2 * pnorm(-abs(tab$meta_z)); tab$meta_FDR <- p.adjust(tab$meta_p, "BH")
tab <- tab[order(tab$disc_p), ]; write.csv(tab, file.path(out, "H1_primary_all_genes.csv"), row.names = FALSE)
cat("Stouffer meta FDR<0.05 (exploratory):", sum(tab$meta_FDR < .05), "\n"); print(head(tab[, c("symbol", "disc_t", "disc_FDR", "rep_t", "rep_p", "meta_FDR")], 15))

# ================= GSVA scores =================
gmt <- strsplit(readLines(gmtf), "\t"); gs <- setNames(lapply(gmt, `[`, -(1:2)), sub("\\|.*", "", sapply(gmt, `[`, 1)))
score <- function(p, s) { E <- P[[p]]$e[genes, s$sample_id]; rownames(E) <- gene_sym
  sc <- gsva(gsvaParam(E, gs, minSize = 10), verbose = FALSE); t(scale(t(sc))) }      # z-scored within platform
cov_gs <- sapply(gs, function(g) sum(g %in% gene_sym)); write.csv(data.frame(set = names(gs), size = lengths(gs), in_universe = cov_gs), file.path(out, "geneset_coverage.csv"), row.names = FALSE)
SC <- list(); for (p in names(P)) { s <- get(cohorts$primary_QC, p); SC[[p]] <- cbind(s, t(score(p, s))) }

# ================= H2 =================
h2sets <- grep("^H2_", names(gs), value = TRUE); comp4 <- c("Macrophage", "AlvFibroblast", "Neutrophil", "Ciliated")
h2 <- data.frame()
for (p in names(P)) for (g in h2sets) { s <- SC[[p]]
  fa <- summary(lm(reformulate(c("DLCO", "FEV1", "sqrtEmph", "Age", "Male", "Current"), g), s))$coef["DLCO", ]
  f0 <- summary(lm(reformulate(c("DLCO", "Age", "Male", "Current"), g), s))$coef["DLCO", ]
  f5 <- summary(lm(reformulate(c("DLCO", "FEV1", "sqrtEmph", "Age", "Male", "Current", comp4), g), s))$coef["DLCO", ]
  h2 <- rbind(h2, data.frame(platform = p, set = g, beta_SD_per10_adj = fa[1] * 10, se = fa[2] * 10, p_adj = fa[4],
    p_one_sided_pos = pt(fa[3], df = nrow(s) - 7, lower.tail = FALSE), beta_per10_unadjFEV1emph = f0[1] * 10, p_unadjFEV1emph = f0[4],
    H5_beta_per10_compadj = f5[1] * 10, H5_p = f5[4])) }
h2$p_adj_holm <- ave(h2$p_adj, h2$platform, FUN = function(x) p.adjust(x, "holm"))
write.csv(h2, file.path(out, "H2_genesets.csv"), row.names = FALSE); print(h2, digits = 3)

# ================= H4 =================
clin <- c("FEV1", "sqrtEmph", "Age", "Male", "Current")
tr <- SC$GPL14550; te <- SC$GPL6480
sig_genes <- order(D$p)[1:100]; w <- sign(D$t[sig_genes])                        # exploratory discovery-derived signature
sigscore <- function(p, s) { E <- P[[p]]$e[genes[sig_genes], s$sample_id]; as.numeric(scale(colMeans(t(scale(t(E))) * w))) }
tr$Signature100 <- sigscore("GPL14550", tr); te$Signature100 <- sigscore("GPL6480", te)
r2oos <- function(f, d) 1 - sum((d$DLCO - predict(f, d))^2) / sum((d$DLCO - mean(d$DLCO))^2)
models <- list(A_clinical = clin, B_clinical_plus_H2sets = c(clin, h2sets), C_clinical_plus_signature100 = c(clin, "Signature100"))
fits <- lapply(models, function(v) lm(reformulate(v, "DLCO"), tr))
bt <- replicate(2000, { i <- sample(nrow(te), replace = TRUE); d <- te[i, ]; sapply(fits, r2oos, d = d) })
h4 <- data.frame(model = names(models), R2_train = sapply(fits, function(f) summary(f)$r.squared),
  R2_test = sapply(fits, r2oos, d = te), R2_test_CI_low = apply(bt, 1, quantile, .025), R2_test_CI_high = apply(bt, 1, quantile, .975))
h4$dR2_test_vs_A <- h4$R2_test - h4$R2_test[1]
dB <- bt[2, ] - bt[1, ]; dC <- bt[3, ] - bt[1, ]
h4$dR2_CI <- c(NA, paste(round(quantile(dB, c(.025, .975)), 3), collapse = " to "), paste(round(quantile(dC, c(.025, .975)), 3), collapse = " to "))
h4$P_dR2_gt0_boot <- c(NA, mean(dB > 0), mean(dC > 0))
write.csv(h4, file.path(out, "H4_prediction.csv"), row.names = FALSE); print(h4, digits = 3)

# ================= H5 (gene level) =================
h5 <- data.frame()
for (p in names(P)) { s <- SC[[p]]; E <- P[[p]]$e[genes, s$sample_id]
  base <- fitD(E, s); adj <- fitD(E, s, update(FORM, as.formula(paste("~ . +", paste(c(h2sets, comp4), collapse = "+")))))
  adj4 <- fitD(E, s, update(FORM, as.formula(paste("~ . +", paste(comp4, collapse = "+")))))
  top <- order(D$p)[1:100]
  for (nm in c("all8_celltypes", "nonalveolar4")) { a <- if (nm == "all8_celltypes") adj else adj4
    h5 <- rbind(h5, data.frame(platform = p, adjustment = nm, FDR05_after = sum(p.adjust(a$p, "BH") < .05),
      tcor_vs_unadjusted_allgenes = cor(base$t, a$t), top100disc_same_direction = mean(sign(a$b[top]) == sign(base$b[top])),
      top100disc_median_beta_ratio = median(a$b[top] / base$b[top]))) } }
write.csv(h5, file.path(out, "H5_composition_adjustment.csv"), row.names = FALSE); print(h5, digits = 3)

# ================= Figures =================
png(file.path(out, "fig5_H2_forest.png"), width = 2400, height = 1300, res = 300)
par(mar = c(4.2, 7.5, 2.2, 1), las = 1, bty = "l")
y <- rev(seq_len(nrow(h2))); col <- ifelse(h2$platform == "GPL14550", "#2a78d6", "#eb6834")
lo <- h2$beta_SD_per10_adj - 1.96 * h2$se; hi <- h2$beta_SD_per10_adj + 1.96 * h2$se
plot(h2$beta_SD_per10_adj, y, xlim = range(lo, hi, 0), pch = 16, col = col, yaxt = "n", ylab = "", xlab = "GSVA score (SD) per +10 %pred DLCO, adjusted")
segments(lo, y, hi, y, col = col, lwd = 2); abline(v = 0, lty = 3, col = "#6b6b68")
axis(2, y, paste(sub("H2_", "", h2$set), ifelse(h2$platform == "GPL14550", "(disc)", "(rep)")), cex.axis = 0.8)
title(main = "H2: pre-specified alveolar cell-type scores (QC cohort)", adj = 0, cex.main = 0.9, font.main = 1)
dev.off()
saveRDS(list(SC = SC, h1 = h1sum, h2 = h2, h3 = h3, h4 = h4, h5 = h5), file.path(out, "results.rds"))

# ---- exploratory: inverse-variance meta of H2 effects across platforms (estimates combined, data not merged) ----
h2m <- do.call(rbind, lapply(split(h2, h2$set), function(d) { w <- 1 / d$se^2; b <- sum(w * d$beta_SD_per10_adj) / sum(w); se <- sqrt(1 / sum(w))
  Q <- sum(w * (d$beta_SD_per10_adj - b)^2)
  data.frame(set = d$set[1], beta_disc = d$beta_SD_per10_adj[1], beta_rep = d$beta_SD_per10_adj[2], meta_beta_per10 = b, meta_se = se,
             meta_p = 2 * pnorm(-abs(b / se)), heterogeneity_Q_p = pchisq(Q, 1, lower.tail = FALSE)) }))
h2m$meta_p_holm <- p.adjust(h2m$meta_p, "holm")
write.csv(h2m, file.path(out, "H2_meta_exploratory.csv"), row.names = FALSE); print(h2m, digits = 3)
