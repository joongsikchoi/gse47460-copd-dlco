# Step 13b: alveolar-capillary fraction ~ DLCO in all 186 (outliers kept) vs QC 177, per platform + IVW meta.
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; ddir <- args[3]; out <- args[4]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
dp <- read.csv(file.path(ddir, "deconv_proportions_186.csv")); dp$AlvCap <- 100 * (dp$AT1 + dp$AT2 + dp$CapA + dp$CapB)
dp$CapA <- 100 * dp$CapA; dp$AT1 <- 100 * dp$AT1
m <- merge(main, dp[, c("sample_id", "AlvCap", "CapA", "AT1")], by = "sample_id")
r <- data.frame()
for (coh in c("all186", "QC177")) for (y in c("AlvCap", "CapA", "AT1")) {
  f <- lapply(c("GPL14550", "GPL6480"), function(p) { s <- m[m$platform == p & (coh == "all186" | !m$qc_outlier), ]
    summary(lm(reformulate(c("DLCO", "FEV1", "sqrtEmph", "Age", "Male", "Current"), y), s))$coef["DLCO", ] })
  b <- sapply(f, `[`, 1) * 10; se <- sapply(f, `[`, 2) * 10; w <- 1 / se^2; mb <- sum(w * b) / sum(w); ms <- sqrt(1 / sum(w))
  r <- rbind(r, data.frame(cohort = coh, y = y, disc_b = b[1], disc_p = f[[1]][4], rep_b = b[2], rep_p = f[[2]][4], meta_b = mb, meta_p = 2 * pnorm(-abs(mb / ms)))) }
write.csv(r, file.path(out, "AlvCap_186_vs_177.csv"), row.names = FALSE); print(r, digits = 3)
