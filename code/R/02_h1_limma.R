# GSE47460 COPD DLCO study - Step 2: H1 gene-wise association (limma) + H3 replication preview
# Model (per platform, never merged):
#   expression ~ DLCO + FEV1 + sqrt(emphysema) + age + sex + current smoking
# Decisions (2026-10-20): main cohort n=186; emphysema sqrt-transformed; smoking = current vs other;
#   GOLD 0 kept in main, excluded in sensitivity.
# Usage: Rscript 02_h1_limma.R <prepared dir> <01_cohort dir> <out dir>
suppressPackageStartupMessages(library(limma))
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; out <- args[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(20261020)

load_platform <- function(p) {
  e <- as.matrix(read.delim(gzfile(file.path(prep, paste0(p, "_expression.tsv.gz"))), row.names = 1, check.names = FALSE))
  a <- read.delim(gzfile(file.path(prep, paste0(p, "_gene_annotation.tsv.gz"))), quote = "", colClasses = "character")
  list(e = e, sym = a$GENE_SYMBOL[match(rownames(e), a$ID)])
}
P <- list(GPL14550 = load_platform("GPL14550"), GPL6480 = load_platform("GPL6480"))
stopifnot(identical(rownames(P$GPL14550$e), rownames(P$GPL6480$e)))

prep_cov <- function(f) {
  d <- read.delim(f, check.names = FALSE)
  d$sqrtEmph <- sqrt(d$Emphysema)
  d$Male <- as.integer(d$Sex == "Male")
  d$Current <- as.integer(d$Smoking == "Current")
  d
}
main <- prep_cov(file.path(cdir, "cohort_main_186.tsv"))
sens188 <- prep_cov(file.path(cdir, "cohort_sensitivity_188.tsv"))

# ---- gene filter (pre-specified, defined on discovery cohort samples only) ----
sym <- P$GPL14550$sym
disc_ids <- main$sample_id[main$platform == "GPL14550"]
mean_disc <- rowMeans(P$GPL14550$e[, disc_ids])
keep <- !is.na(sym) & sym != "" & mean_disc > quantile(mean_disc, 0.25)
# one row per symbol (highest discovery mean)
ord <- order(-mean_disc); kk <- ord[keep[ord]]; kk <- kk[!duplicated(sym[kk])]
genes <- rownames(P$GPL14550$e)[kk]
cat("Genes analysed:", length(genes), "\n")

fit_platform <- function(d, p, form = ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current) {
  s <- d[d$platform == p, ]
  E <- P[[p]]$e[genes, s$sample_id]
  X <- model.matrix(form, s)
  fit <- eBayes(lmFit(E, X), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "DLCO", number = Inf, sort.by = "none")
  tt$symbol <- P[[p]]$sym[match(rownames(tt), rownames(P[[p]]$e))]
  list(tt = tt, fit = fit, E = E, X = X, n = nrow(s))
}
vif <- function(s) { X <- s[, c("DLCO", "FEV1", "sqrtEmph", "Age", "Male", "Current")]
  sapply(names(X), function(v) 1 / (1 - summary(lm(reformulate(setdiff(names(X), v), v), X))$r.squared)) }

res <- list(); log <- c()
analyses <- list(
  main          = list(d = main,  form = ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current),
  no_age_sex    = list(d = main,  form = ~ DLCO + FEV1 + sqrtEmph + Current),
  excl_GOLD0    = list(d = main[is.na(main$GOLD) | main$GOLD != "At Risk", ], form = ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current),
  cohort188     = list(d = sens188, form = ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current),
  no_FEV1_emph  = list(d = main,  form = ~ DLCO + Age + Male + Current)   # reference: unadjusted for lung structure
)
summ <- data.frame()
for (an in names(analyses)) for (p in names(P)) {
  r <- fit_platform(analyses[[an]]$d, p, analyses[[an]]$form)
  res[[an]][[p]] <- r
  tt <- r$tt
  summ <- rbind(summ, data.frame(analysis = an, platform = p, n = r$n,
    FDR05 = sum(tt$adj.P.Val < 0.05), FDR10 = sum(tt$adj.P.Val < 0.10), FDR25 = sum(tt$adj.P.Val < 0.25),
    nominal_p05 = sum(tt$P.Value < 0.05), expected_p05 = round(0.05 * nrow(tt)),
    pi0 = round(propTrueNull(tt$P.Value), 3)))
}
write.csv(summ, file.path(out, "h1_summary_counts.csv"), row.names = FALSE)
print(summ)

for (p in names(P)) {
  v <- vif(main[main$platform == p, ]); cat("VIF", p, round(v, 2), "\n")
  log <- c(log, paste("VIF", p, paste(names(v), round(v, 2), collapse = "; ")))
}

# ---- full results tables ----
D <- res$main$GPL14550$tt; R <- res$main$GPL6480$tt
full <- data.frame(probe = rownames(D), symbol = D$symbol,
  disc_logFC_per10 = D$logFC * 10, disc_t = D$t, disc_p = D$P.Value, disc_FDR = D$adj.P.Val, disc_AveExpr = D$AveExpr,
  rep_logFC_per10 = R$logFC * 10, rep_t = R$t, rep_p = R$P.Value, rep_FDR = R$adj.P.Val)
full <- full[order(full$disc_p), ]
write.csv(full, file.path(out, "h1_main_results_all_genes.csv"), row.names = FALSE)

# ---- H3 preview: replication of discovery signals (Freedman-Lane permutation on DLCO in replication) ----
rep_perm <- function(sel, B = 1000) {
  r <- res$main$GPL6480; X <- r$X
  Z <- X[, colnames(X) != "DLCO"]
  dl_hat <- fitted(lm(X[, "DLCO"] ~ Z - 1)); dl_res <- X[, "DLCO"] - dl_hat
  dsign <- sign(D[sel, "logFC"])
  stat <- function(tvec) c(conc = mean(sign(tvec[sel]) == dsign), cor = cor(D[sel, "t"], tvec[sel]))
  obs <- stat(setNames(R$t, rownames(R))); names(obs) <- c("conc", "cor")
  E <- r$E[sel, , drop = FALSE]
  null <- t(replicate(B, { Xp <- X; Xp[, "DLCO"] <- dl_hat + sample(dl_res)
    f <- eBayes(lmFit(E, Xp), trend = TRUE); stat(setNames(f$t[, "DLCO"], rownames(E))[sel]) }))
  c(n_genes = length(sel), concordance = obs[["conc"]], p_conc = (1 + sum(null[, 1] >= obs[["conc"]])) / (B + 1),
    t_correlation = obs[["cor"]], p_cor = (1 + sum(null[, 2] >= obs[["cor"]])) / (B + 1))
}
sets <- list(FDR05 = rownames(D)[D$adj.P.Val < 0.05], FDR10 = rownames(D)[D$adj.P.Val < 0.10],
             top100 = rownames(D)[order(D$P.Value)][1:100], top500 = rownames(D)[order(D$P.Value)][1:500],
             all_genes = rownames(D))
h3 <- do.call(rbind, lapply(names(sets), function(s) {
  sel <- sets[[s]]; if (length(sel) < 5) return(data.frame(set = s, n_genes = length(sel)))
  data.frame(set = s, t(rep_perm(sel, B = if (s == "all_genes") 200 else 1000)))
}))
write.csv(h3, file.path(out, "h3_replication_preview.csv"), row.names = FALSE)
print(h3)

# replicated genes: discovery FDR<0.05 (or FDR<0.10) + same direction + replication p<0.05 one-sided
full$replicated <- with(full, disc_FDR < 0.10 & sign(disc_t) == sign(rep_t) & rep_p / 2 < 0.05)
write.csv(full[full$disc_FDR < 0.10, ], file.path(out, "h1_discovery_FDR10_genes.csv"), row.names = FALSE)

# ---- robustness of discovery signal across sensitivity analyses ----
rob <- sapply(c("main", "no_age_sex", "excl_GOLD0", "cohort188", "no_FEV1_emph"), function(an) res[[an]]$GPL14550$tt[rownames(D), "t"])
write.csv(round(cor(rob), 3), file.path(out, "h1_sensitivity_t_correlations.csv"))

# ---- figures ----
png(file.path(out, "fig2_h1_volcano_qq.png"), width = 3600, height = 1500, res = 300)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.5, 1), las = 1, bty = "l", col.axis = "#4a4a48", col.lab = "#1f1f1e")
for (p in names(P)) {
  tt <- res$main[[p]]$tt; sig <- tt$adj.P.Val < 0.05
  plot(tt$logFC * 10, -log10(tt$P.Value), pch = 16, cex = 0.4, col = ifelse(sig, ifelse(tt$logFC > 0, "#2a78d6", "#eb6834"), "#bdbdb8"),
       xlab = "log2 FC per +10 %pred DLCO", ylab = "-log10 p", main = "")
  title(main = sprintf("%s (n=%d): %d genes FDR<0.05", p, res$main[[p]]$n, sum(sig)), adj = 0, cex.main = 1, font.main = 1)
  top <- head(order(tt$P.Value), 8)
  text(tt$logFC[top] * 10, -log10(tt$P.Value[top]), tt$symbol[top], pos = 3, cex = 0.6, col = "#1f1f1e")
}
pv <- res$main$GPL14550$tt$P.Value; pr <- res$main$GPL6480$tt$P.Value; ex <- -log10(ppoints(length(pv)))
plot(ex, sort(-log10(pv), decreasing = TRUE), pch = 16, cex = 0.4, col = "#2a78d6", xlab = "Expected -log10 p", ylab = "Observed -log10 p",
     ylim = range(0, -log10(c(pv, pr))))
points(ex, sort(-log10(pr), decreasing = TRUE), pch = 16, cex = 0.4, col = "#eb6834"); abline(0, 1, col = "#6b6b68", lty = 2)
legend("topleft", c("GPL14550 (discovery)", "GPL6480 (replication)"), col = c("#2a78d6", "#eb6834"), pch = 16, bty = "n", cex = 0.8)
title(main = "QQ plot, DLCO coefficient", adj = 0, cex.main = 1, font.main = 1)
dev.off()

png(file.path(out, "fig3_h3_t_concordance.png"), width = 1800, height = 1650, res = 300)
par(mar = c(4.2, 4.2, 2.5, 1), las = 1, bty = "l")
top500 <- sets$top500
plot(D$t, R$t, pch = 16, cex = 0.3, col = "#d9d9d5", xlab = "Discovery t (GPL14550)", ylab = "Replication t (GPL6480)")
points(D[top500, "t"], R[top500, "t"], pch = 16, cex = 0.5, col = "#2a78d6"); abline(h = 0, v = 0, col = "#6b6b68", lty = 3)
title(main = sprintf("DLCO t-statistics; top 500 discovery genes: r=%.2f", cor(D[top500, "t"], R[top500, "t"])), adj = 0, cex.main = 0.9, font.main = 1)
dev.off()
writeLines(c(paste("Genes analysed:", length(genes)), log), file.path(out, "h1_log.txt"))

# ---- contrast: cross-platform agreement for the model WITHOUT FEV1/emphysema adjustment ----
Du <- res$no_FEV1_emph$GPL14550$tt; Ru <- res$no_FEV1_emph$GPL6480$tt
topu <- order(Du$P.Value)[1:500]
contrast <- data.frame(model = c("adjusted (main)", "unadjusted for FEV1/emphysema"),
  t_cor_all = c(cor(D$t, R$t), cor(Du$t, Ru$t)),
  t_cor_top500 = c(cor(D$t[order(D$P.Value)[1:500]], R$t[order(D$P.Value)[1:500]]), cor(Du$t[topu], Ru$t[topu])),
  concordance_top500 = c(mean(sign(D$t[order(D$P.Value)[1:500]]) == sign(R$t[order(D$P.Value)[1:500]])), mean(sign(Du$t[topu]) == sign(Ru$t[topu]))))
write.csv(contrast, file.path(out, "h3_adjusted_vs_unadjusted_agreement.csv"), row.names = FALSE); print(contrast)
