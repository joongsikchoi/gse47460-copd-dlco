# Step 13e: robust regression (Huber M-estimation, MASS::rlm) of estimated cell fractions on DLCO
# in all 186 patients (outliers kept) and in the main analysis (177), per platform + IVW meta-analysis.
# Requested in the 2026-11-10 review. Same model and inputs as 13b_alvcap_186.R.
# Usage: Rscript 13e_rlm_186.R <01_cohort dir> <10_deconvolution dir> <out dir>
suppressPackageStartupMessages(library(MASS))
args <- commandArgs(TRUE); cdir <- args[1]; ddir <- args[2]; out <- args[3]; dir.create(out, FALSE, TRUE)
d <- read.delim(file.path(cdir, "cohort_main_186.tsv"), check.names = FALSE)
d$sqrtEmph <- sqrt(d$Emphysema); d$Male <- as.integer(d$Sex == "Male"); d$Current <- as.integer(d$Smoking == "Current")
d$qc_outlier <- d$sample_id %in% readLines(file.path(cdir, "qc_outliers_9.txt"))
dp <- read.csv(file.path(ddir, "deconv_proportions_186.csv"))
dp$AlvCap <- 100 * (dp$AT1 + dp$AT2 + dp$CapA + dp$CapB); dp$CapA <- 100 * dp$CapA; dp$AT1 <- 100 * dp$AT1
m <- merge(d, dp[, c("sample_id", "AlvCap", "CapA", "AT1")], by = "sample_id")
stopifnot(nrow(m) == 186, sum(m$qc_outlier) == 9)
X <- c("DLCO", "FEV1", "sqrtEmph", "Age", "Male", "Current")
r <- data.frame(); wt <- data.frame()
for (coh in c("all186", "QC177")) for (y in c("AlvCap", "CapA", "AT1")) for (method in c("lm", "rlm_huber")) {
  f <- lapply(c("GPL14550", "GPL6480"), function(p) {
    s <- m[m$platform == p & (coh == "all186" | !m$qc_outlier), ]
    if (method == "lm") { co <- summary(lm(reformulate(X, y), s))$coef["DLCO", ]; c(b = co[[1]], se = co[[2]], p = co[[4]], n = nrow(s), conv = 1) }
    else { fit <- rlm(reformulate(X, y), s, psi = psi.huber, maxit = 200); co <- summary(fit)$coef["DLCO", ]
      if (coh == "all186" && y == "AlvCap") wt <<- rbind(wt, data.frame(sample_id = s$sample_id, platform = p, qc_outlier = s$qc_outlier, AlvCap = s$AlvCap, huber_weight = fit$w))
      c(b = co[[1]], se = co[[2]], p = 2 * pt(-abs(co[[3]]), nrow(s) - length(X) - 1), n = nrow(s), conv = as.numeric(fit$converged)) } })
  b <- sapply(f, `[`, "b") * 10; se <- sapply(f, `[`, "se") * 10; w <- 1 / se^2; mb <- sum(w * b) / sum(w); ms <- sqrt(1 / sum(w))
  r <- rbind(r, data.frame(cohort = coh, y = y, method = method, disc_n = f[[1]][["n"]], rep_n = f[[2]][["n"]],
    disc_b = b[[1]], disc_p = f[[1]][["p"]], rep_b = b[[2]], rep_p = f[[2]][["p"]],
    meta_b = mb, meta_lo = mb - 1.96 * ms, meta_hi = mb + 1.96 * ms, meta_p = 2 * pnorm(-abs(mb / ms)),
    converged = f[[1]][["conv"]] == 1 && f[[2]][["conv"]] == 1)) }
write.csv(r, file.path(out, "I_rlm_186_vs_177.csv"), row.names = FALSE)
write.csv(wt, file.path(out, "I_rlm_huber_weights_AlvCap_186.csv"), row.names = FALSE)
print(r, digits = 3)
cat("\nHuber weights, AlvCap model, all 186 (outliers):\n"); print(wt[wt$qc_outlier, ], digits = 3)
cat("\nNumber of arrays with weight < 1:", sum(wt$huber_weight < 1), "of", nrow(wt), "\n")
cat("Lowest weights:\n"); print(head(wt[order(wt$huber_weight), ], 8), digits = 3)
