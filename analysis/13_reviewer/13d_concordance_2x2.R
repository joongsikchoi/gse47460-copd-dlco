# Step 13d: isolate which platform's outlier arrays drive low concordance; share of the replicated signal that tracks AlvCap.
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; ddir <- args[3]; out <- args[4]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
dp <- read.csv(file.path(ddir, "deconv_proportions_186.csv")); dp$AlvCap <- 100 * (dp$AT1 + dp$AT2 + dp$CapA + dp$CapB)
main <- merge(main, dp[, c("sample_id", "AlvCap")], by = "sample_id", sort = FALSE)
FA <- update(FORM, ~ . + AlvCap)
print(table(main$platform, main$qc_outlier))
fp <- function(p, qc, form) { s <- main[main$platform == p & (if (qc) !main$qc_outlier else TRUE), ]; fitD(P[[p]]$e[genes, s$sample_id], s, form) }
conc <- function(D, R, k = 500) { s <- order(D$p)[1:k]; mean(sign(D$t[s]) == sign(R$t[s])) }
r <- data.frame()
for (form in c("base", "AlvCap")) { F <- if (form == "base") FORM else FA
  D <- list(QC = fp("GPL14550", TRUE, F), all = fp("GPL14550", FALSE, F)); R <- list(QC = fp("GPL6480", TRUE, F), all = fp("GPL6480", FALSE, F))
  for (d in names(D)) for (q in names(R)) r <- rbind(r, data.frame(model = form, discovery = d, replication = q, conc_top500 = conc(D[[d]], R[[q]]), conc_top100 = conc(D[[d]], R[[q]], 100))) }
print(r); write.csv(r, file.path(out, "H2_concordance_2x2.csv"), row.names = FALSE)
# how much of the QC177 base top-500 replication rides on composition
D <- fp("GPL14550", TRUE, FORM); R <- fp("GPL6480", TRUE, FORM); sel <- order(D$p)[1:500]
sD <- main[main$platform == "GPL14550" & !main$qc_outlier, ]; ra <- apply(P$GPL14550$e[genes[sel], sD$sample_id], 1, cor, y = sD$AlvCap)
cat("QC177 base top500: sign(t) agrees with sign(cor AlvCap):", round(mean(sign(D$t[sel]) == sign(ra)), 3), " median |r|:", round(median(abs(ra)), 3), "\n")
hi <- abs(ra) >= median(abs(ra))
cat("concordance in genes with |r_AlvCap| above median:", round(mean(sign(D$t[sel][hi]) == sign(R$t[sel][hi])), 3), " below:", round(mean(sign(D$t[sel][!hi]) == sign(R$t[sel][!hi])), 3), "\n")
