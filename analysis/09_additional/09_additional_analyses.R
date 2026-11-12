# Step 9: additional analyses for the manuscript (QC cohort 118 / 59), all pre-planned in the protocol or requested as robustness checks.
#  A. model variants: without age/sex (protocol: %pred already reflects age/sex), + upper-lobe indicator
#  B. H2 alternative scoring (mean of within-platform gene z-scores instead of GSVA)
#  C. H2 meta leave-one-out influence
#  D. auxiliary outcome: lowest tertile of residual DLCO (protocol 2-B) - group comparison + discovery-trained classifier AUC in replication
#  E. H2 scores vs FEV1 / emphysema / raw DLCO (interpretation)
#  F. competitive gene-set test (limma camera) on GO biological processes, per platform, with cross-platform agreement
# Usage: Rscript 09_additional_analyses.R <prepared dir> <01_cohort dir> <gmt> <04_main out dir> <out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; gmtf <- args[3]; o4 <- args[4]; out <- args[5]; dir.create(out, FALSE, TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
suppressPackageStartupMessages({ library(GSVA); library(pROC) }); set.seed(20261029)
R4 <- readRDS(file.path(o4, "results.rds")); SC <- R4$SC
gmt <- strsplit(readLines(gmtf), "\t"); gs <- setNames(lapply(gmt, `[`, -(1:2)), sub("\\|.*", "", sapply(gmt, `[`, 1)))
h2sets <- grep("^H2_", names(gs), value = TRUE); clin <- c("FEV1", "sqrtEmph", "Age", "Male", "Current")
ivw <- function(b, se) { w <- 1 / se^2; m <- sum(w * b) / sum(w); s <- sqrt(1 / sum(w)); c(meta = m, meta_se = s, meta_p = 2 * pnorm(-abs(m / s)), Q_p = pchisq(sum(w * (b - m)^2), 1, lower.tail = FALSE)) }
h2fit <- function(s, covs, score_cols = h2sets) t(sapply(score_cols, function(g) summary(lm(reformulate(c("DLCO", covs), g), s))$coef["DLCO", c(1, 2, 4)] * c(10, 10, 1)))
conc100 <- function(D, R) { top <- order(D$p)[1:100]; mean(sign(D$t[top]) == sign(R$t[top])) }

# ---------- A. model variants ----------
variants <- list(main = clin, no_age_sex = c("FEV1", "sqrtEmph", "Current"), plus_upper_lobe = c(clin, "UpperLobe"))
A <- data.frame()
for (v in names(variants)) {
  form <- reformulate(c("DLCO", variants[[v]])); sD <- SC$GPL14550; sR <- SC$GPL6480
  D <- fitD(P$GPL14550$e[genes, sD$sample_id], sD, form); R <- fitD(P$GPL6480$e[genes, sR$sample_id], sR, form)
  hD <- h2fit(sD, variants[[v]]); hR <- h2fit(sR, variants[[v]])
  row <- data.frame(model = v, disc_FDR05 = sum(p.adjust(D$p, "BH") < .05), rep_FDR05 = sum(p.adjust(R$p, "BH") < .05), conc_top100 = conc100(D, R))
  for (g in h2sets) { m <- ivw(c(hD[g, 1], hR[g, 1]), c(hD[g, 2], hR[g, 2]))
    row[[paste0(sub("H2_", "", g), "_disc_p")]] <- hD[g, 3]; row[[paste0(sub("H2_", "", g), "_meta_beta")]] <- m["meta"]; row[[paste0(sub("H2_", "", g), "_meta_p")]] <- m["meta_p"] }
  A <- rbind(A, row)
}
write.csv(A, file.path(out, "A_model_variants.csv"), row.names = FALSE); print(A, digits = 3)

# ---------- B. alternative score: mean z ----------
B <- data.frame()
for (p in names(P)) { s <- SC[[p]]; E <- P[[p]]$e[genes, s$sample_id]; rownames(E) <- gene_sym; Z <- t(scale(t(E)))
  for (g in h2sets) { gg <- intersect(gs[[g]], rownames(Z)); s[[paste0("z_", g)]] <- as.numeric(scale(colMeans(Z[gg, ]))) }
  SC[[p]] <- s
  hz <- h2fit(s, clin, paste0("z_", h2sets))
  B <- rbind(B, data.frame(platform = p, set = sub("H2_", "", h2sets), beta = hz[, 1], se = hz[, 2], p = hz[, 3],
                           cor_with_gsva = sapply(h2sets, function(g) cor(s[[g]], s[[paste0("z_", g)]])))) }
Bm <- do.call(rbind, lapply(split(B, B$set), function(d) data.frame(set = d$set[1], t(ivw(d$beta, d$se)))))
write.csv(B, file.path(out, "B_meanz_scores.csv"), row.names = FALSE); write.csv(Bm, file.path(out, "B_meanz_meta.csv"), row.names = FALSE)
print(B, digits = 3); print(Bm, digits = 3)

# ---------- C. leave-one-out influence on H2 meta ----------
C <- data.frame()
for (g in c("H2_AT1", "H2_gCap", "H2_aCap")) for (p in names(P)) { s <- SC[[p]]; other <- h2fit(SC[[setdiff(names(P), p)]], clin, g)
  for (i in seq_len(nrow(s))) { h <- h2fit(s[-i, ], clin, g); m <- ivw(c(h[1, 1], other[1, 1]), c(h[1, 2], other[1, 2]))
    C <- rbind(C, data.frame(set = sub("H2_", "", g), dropped_platform = p, dropped = s$sample_id[i], meta_beta = m["meta"], meta_p = m["meta_p"])) } }
Cs <- do.call(rbind, lapply(split(C, C$set), function(d) data.frame(set = d$set[1], n_loo = nrow(d), meta_p_max = max(d$meta_p),
  dropped_at_max = d$dropped[which.max(d$meta_p)], meta_beta_min = min(d$meta_beta), meta_beta_max = max(d$meta_beta), frac_p_below_05 = mean(d$meta_p < .05))))
write.csv(Cs, file.path(out, "C_loo_H2_meta.csv"), row.names = FALSE); print(Cs, digits = 3)

# ---------- D. auxiliary binary outcome: lowest tertile of residual DLCO ----------
for (p in names(P)) { s <- SC[[p]]; s$residDLCO <- resid(lm(reformulate(clin, "DLCO"), s))
  s$lowT <- as.integer(s$residDLCO <= quantile(s$residDLCO, 1/3)); SC[[p]] <- s }
Dg <- data.frame()
for (p in names(P)) for (g in h2sets) { s <- SC[[p]]; w <- wilcox.test(s[[g]][s$lowT == 1], s[[g]][s$lowT == 0])
  Dg <- rbind(Dg, data.frame(platform = p, set = sub("H2_", "", g), n_low = sum(s$lowT), median_low = median(s[[g]][s$lowT == 1]),
    median_rest = median(s[[g]][s$lowT == 0]), diff_mean_SD = mean(s[[g]][s$lowT == 1]) - mean(s[[g]][s$lowT == 0]), wilcox_p = w$p.value)) }
tr <- SC$GPL14550; te <- SC$GPL6480
fit <- glm(reformulate(h2sets, "lowT"), binomial, tr); pr <- predict(fit, te, type = "response")
roc_te <- roc(te$lowT, pr, quiet = TRUE, direction = "<"); ci <- ci.auc(roc_te, method = "bootstrap", boot.n = 2000)
roc_tr <- roc(tr$lowT, fitted(fit), quiet = TRUE, direction = "<")
fit1 <- glm(lowT ~ H2_AT1, binomial, tr); roc1 <- roc(te$lowT, predict(fit1, te, type = "response"), quiet = TRUE, direction = "<"); ci1 <- ci.auc(roc1, method = "bootstrap", boot.n = 2000)
Dc <- data.frame(model = c("4 H2 scores", "AT1 score only"), AUC_train = c(auc(roc_tr), auc(roc(tr$lowT, fitted(fit1), quiet = TRUE, direction = "<"))),
  AUC_replication = c(auc(roc_te), auc(roc1)), CI_low = c(ci[1], ci1[1]), CI_high = c(ci[3], ci1[3]))
write.csv(Dg, file.path(out, "D_lowtertile_groups.csv"), row.names = FALSE); write.csv(Dc, file.path(out, "D_lowtertile_classifier.csv"), row.names = FALSE)
print(Dg, digits = 3); print(Dc, digits = 3)
sav <- do.call(rbind, lapply(names(P), function(p) data.frame(platform = p, SC[[p]][, c("sample_id", "DLCO", "residDLCO", "lowT", h2sets)])))
write.csv(sav, file.path(out, "D_scores_by_tertile.csv"), row.names = FALSE)

# ---------- E. H2 scores vs clinical measures (Spearman) ----------
E5 <- do.call(rbind, lapply(names(P), function(p) { s <- SC[[p]]; do.call(rbind, lapply(h2sets, function(g) data.frame(platform = p, set = sub("H2_", "", g),
  rho_DLCO = cor(s[[g]], s$DLCO, method = "s"), rho_FEV1 = cor(s[[g]], s$FEV1, method = "s"), rho_emph = cor(s[[g]], s$Emphysema, method = "s"),
  rho_residDLCO = cor(s[[g]], s$residDLCO, method = "s")))) }))
write.csv(E5, file.path(out, "E_scores_vs_clinical.csv"), row.names = FALSE); print(E5, digits = 2)

# ---------- F. camera on GO biological process ----------
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(GO.db) })
go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(gene_sym), keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL")))
go <- go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", c("SYMBOL", "GOALL")]; go <- unique(go)
idx_all <- split(match(go$SYMBOL, gene_sym), go$GOALL); idx_all <- lapply(idx_all, function(i) unique(i[!is.na(i)]))
idx_all <- idx_all[lengths(idx_all) >= 15 & lengths(idx_all) <= 300]
cam <- list()
for (p in names(P)) { s <- SC[[p]]; X <- model.matrix(FORM, s); E <- P[[p]]$e[genes, s$sample_id]
  cam[[p]] <- camera(E, idx_all, X, contrast = "DLCO", inter.gene.cor = 0.01, sort = FALSE) }
F6 <- data.frame(GO = names(idx_all), term = suppressMessages(AnnotationDbi::select(GO.db, names(idx_all), "TERM", "GOID"))$TERM, size = lengths(idx_all),
  disc_dir = cam$GPL14550$Direction, disc_p = cam$GPL14550$PValue, disc_FDR = cam$GPL14550$FDR,
  rep_dir = cam$GPL6480$Direction, rep_p = cam$GPL6480$PValue, rep_FDR = cam$GPL6480$FDR)
F6$replicated <- F6$disc_FDR < .05 & F6$rep_p < .05 & F6$disc_dir == F6$rep_dir
F6 <- F6[order(F6$disc_p), ]
write.csv(F6, file.path(out, "F_camera_GO_BP.csv"), row.names = FALSE)
topD <- F6$disc_FDR < .05
cat("GO BP sets:", nrow(F6), " disc FDR<.05:", sum(topD), " (Up", sum(topD & F6$disc_dir == "Up"), ") rep FDR<.05:", sum(F6$rep_FDR < .05),
    " replicated (disc FDR<.05 & rep p<.05 same dir):", sum(F6$replicated), "\n")
cat("direction agreement among disc FDR<.05:", mean(F6$disc_dir[topD] == F6$rep_dir[topD]), "\n")
print(head(F6[F6$replicated, c("term", "size", "disc_dir", "disc_FDR", "rep_p")], 30))
# camera on the Travaglini cell-type sets too (competitive counterpart of H2)
ct <- lapply(gs, function(g) which(gene_sym %in% g)); ctab <- do.call(rbind, lapply(names(P), function(p) { s <- SC[[p]]
  r <- camera(P[[p]]$e[genes, s$sample_id], ct, model.matrix(FORM, s), contrast = "DLCO", inter.gene.cor = 0.01, sort = FALSE)
  data.frame(platform = p, set = rownames(r), r) }))
write.csv(ctab, file.path(out, "F_camera_celltype_sets.csv"), row.names = FALSE); print(ctab, digits = 3)
saveRDS(SC, file.path(out, "SC_with_tertile.rds"))
