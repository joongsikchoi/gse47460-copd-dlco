# Step 8a: save split-half draws (direction concordance of top-500 genes) for Figure 3, before and after array QC.
# Same procedure as 03_replication_diagnostics.R (200 random splits, test half size = replication n).
# Usage: Rscript 08a_split_half_draws.R <prepared dir> <01_cohort dir> <out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; out <- args[3]; dir.create(out, FALSE, TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
set.seed(20261020)
conc <- function(a, b, k = 500) { top <- order(a$p)[1:k]; c(conc = mean(sign(a$t[top]) == sign(b$t[top])), tcor = cor(a$t[top], b$t[top])) }
res <- list(); cross <- list()
for (q in c("preQC_186", "QC_177")) {
  keep <- if (q == "QC_177") !main$qc_outlier else rep(TRUE, nrow(main))
  sD <- main[main$platform == "GPL14550" & keep, ]; sR <- main[main$platform == "GPL6480" & keep, ]
  ED <- P$GPL14550$e[genes, sD$sample_id]; ER <- P$GPL6480$e[genes, sR$sample_id]
  cross[[q]] <- data.frame(cohort = q, n_disc = nrow(sD), n_rep = nrow(sR), t(conc(fitD(ED, sD), fitD(ER, sR))))
  res[[q]] <- data.frame(cohort = q, draw = 1:200, t(replicate(200, { idx <- sample(nrow(sD), nrow(sR)); a <- sD[-idx, ]; b <- sD[idx, ]
    conc(fitD(ED[, a$sample_id], a), fitD(ED[, b$sample_id], b)) })))
}
d <- do.call(rbind, res); x <- do.call(rbind, cross)
write.csv(d, file.path(out, "split_half_draws.csv"), row.names = FALSE); write.csv(x, file.path(out, "cross_platform_top500.csv"), row.names = FALSE)
print(x); print(aggregate(cbind(conc, tcor) ~ cohort, d, function(v) round(quantile(v, c(.05, .5, .95)), 3)))
