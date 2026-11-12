# Step 14b: per-platform check of GPIHBP1 and HOPX after own-fraction adjustment (QC177).
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; ddir <- args[3]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
dp <- read.csv(file.path(ddir, "deconv_proportions_186.csv")); main <- merge(main, dp[, c("sample_id", "AT1", "CapB")], by = "sample_id")
for (x in list(c("GPIHBP1", "CapB"), c("HOPX", "AT1"))) for (p in c("GPL14550", "GPL6480")) {
  s <- main[main$platform == p & !main$qc_outlier, ]; s$y <- P[[p]]$e[genes[match(x[1], gene_sym)], s$sample_id]; s$FRAC <- 100 * s[[x[2]]]
  a <- summary(lm(y ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current, s))$coef["DLCO", ]; b <- summary(lm(y ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current + FRAC, s))$coef["DLCO", ]
  cat(x[1], p, "base b/10:", round(10 * a[1], 3), "p", signif(a[4], 2), " | frac-adj b/10:", round(10 * b[1], 3), "p", signif(b[4], 2), " r(y,frac)", round(cor(s$y, s$FRAC), 2), "\n") }
