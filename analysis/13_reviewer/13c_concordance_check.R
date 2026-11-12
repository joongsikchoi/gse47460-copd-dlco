# Step 13c: why top-500 cross-platform concordance is 32.4% in the 186-patient analysis adjusted for the alveolar-capillary fraction.
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; ddir <- args[3]; out <- args[4]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
set.seed(20261020)
dp <- read.csv(file.path(ddir, "deconv_proportions_186.csv")); dp$AlvCap <- 100 * (dp$AT1 + dp$AT2 + dp$CapA + dp$CapB)
main <- merge(main, dp[, c("sample_id", "AlvCap")], by = "sample_id", sort = FALSE)
FA <- update(FORM, ~ . + AlvCap)
fit2 <- function(filt, form) lapply(c("GPL14550", "GPL6480"), function(p) { s <- main[main$platform == p & filt, ]; fitD(P[[p]]$e[genes, s$sample_id], s, form) })
conc <- function(D, R, k) { s <- order(D$p)[1:k]; mean(sign(D$t[s]) == sign(R$t[s])) }
res <- data.frame(); keep <- list()
for (coh in c("all186", "QC177")) for (adj in c("base", "AlvCap")) {
  filt <- if (coh == "all186") rep(TRUE, nrow(main)) else !main$qc_outlier
  f <- fit2(filt, if (adj == "base") FORM else FA); D <- f[[1]]; R <- f[[2]]; keep[[paste(coh, adj)]] <- f
  # permutation null for top-500 concordance: shuffle DLCO residual on the replication platform (Freedman-Lane), B = 200
  s <- main[main$platform == "GPL6480" & filt, ]; X <- model.matrix(if (adj == "base") FORM else FA, s)
  Z <- X[, colnames(X) != "DLCO"]; h <- fitted(lm(X[, "DLCO"] ~ Z - 1)); e <- X[, "DLCO"] - h
  sel <- order(D$p)[1:500]; ER <- P$GPL6480$e[genes[sel], s$sample_id]
  nul <- replicate(200, { s2 <- s; s2$DLCO <- h + sample(e); mean(sign(fitD(ER, s2, if (adj == "base") FORM else FA)$t) == sign(D$t[sel])) })
  res <- rbind(res, data.frame(cohort = coh, adjustment = adj, conc_top100 = conc(D, R, 100), conc_top500 = conc(D, R, 500),
    conc_all_genes = mean(sign(D$t) == sign(R$t)), tcor_all = cor(D$t, R$t), tcor_top500 = cor(D$t[sel], R$t[sel]),
    perm_null_median = median(nul), perm_null_q025 = quantile(nul, .025), perm_null_q975 = quantile(nul, .975),
    perm_p_low = (1 + sum(nul <= conc(D, R, 500))) / 201,
    frac_top500_negative_disc = mean(D$t[sel] < 0)))
}
print(res, digits = 3); write.csv(res, file.path(out, "H_concordance_check.csv"), row.names = FALSE)
# relation of the adjusted discovery top-500 to the alveolar-capillary fraction and to the outliers
D <- keep[["all186 AlvCap"]][[1]]; R <- keep[["all186 AlvCap"]][[2]]; sel <- order(D$p)[1:500]
sD <- main[main$platform == "GPL14550", ]; ED <- P$GPL14550$e[genes[sel], sD$sample_id]
r_alv <- apply(ED, 1, function(x) cor(x, sD$AlvCap))
Dq <- keep[["QC177 AlvCap"]][[1]]
cat("top500 (186+AlvCap): median |cor with AlvCap|", round(median(abs(r_alv)), 3), "\n")
cat("sign(DLCO t) vs sign(cor with AlvCap): agreement", round(mean(sign(D$t[sel]) == sign(r_alv)), 3), "\n")
cat("overlap of top500 with top500 of the QC177+AlvCap analysis:", length(intersect(sel, order(Dq$p)[1:500])), "\n")
cat("same sign as QC177+AlvCap for these genes:", round(mean(sign(D$t[sel]) == sign(Dq$t[sel])), 3), "\n")
# drop the two dominant arrays only, keep AlvCap adjustment
two <- c("GSM1150187", "GSM1150251"); f2 <- fit2(!(main$sample_id %in% two), FA)
cat("186 minus 2 dominant arrays + AlvCap: top500 concordance", round(conc(f2[[1]], f2[[2]], 500), 3), " FDR05 disc", sum(p.adjust(f2[[1]]$p, "BH") < .05), "\n")
