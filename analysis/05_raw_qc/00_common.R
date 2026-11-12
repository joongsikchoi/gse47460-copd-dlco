# Shared loaders for GSE47460 COPD DLCO analyses. Source after setting `prep` (prepared dir) and `cdir` (01_cohort dir).
suppressPackageStartupMessages(library(limma))
load_platform <- function(p) {
  e <- as.matrix(read.delim(gzfile(file.path(prep, paste0(p, "_expression.tsv.gz"))), row.names = 1, check.names = FALSE))
  a <- read.delim(gzfile(file.path(prep, paste0(p, "_gene_annotation.tsv.gz"))), quote = "", colClasses = "character")
  list(e = e, sym = a$GENE_SYMBOL[match(rownames(e), a$ID)])
}
prep_cov <- function(f) {
  d <- read.delim(f, check.names = FALSE)
  d$sqrtEmph <- sqrt(d$Emphysema); d$Male <- as.integer(d$Sex == "Male"); d$Current <- as.integer(d$Smoking == "Current")
  d$lobe <- sub("^LT\\d{6}(..)_.*", "\\1", d$title); d$UpperLobe <- as.integer(substr(d$lobe, 2, 2) == "U")
  d
}
P <- list(GPL14550 = load_platform("GPL14550"), GPL6480 = load_platform("GPL6480"))
main <- prep_cov(file.path(cdir, "cohort_main_186.tsv"))
# sex-discordant samples (RPS4Y1 log2 > 9 vs recorded sex), identified in step 2 QC
sexcheck <- function(p) { e <- P[[p]]$e; v <- e[which(P[[p]]$sym == "RPS4Y1")[1], ]; s <- main[main$platform == p, ]
  s$sample_id[(v[s$sample_id] > 9) != (s$Male == 1)] }
main$sex_discordant <- main$sample_id %in% c(sexcheck("GPL14550"), sexcheck("GPL6480"))
# gene universe (pre-specified in step 2): discovery-mean > 25th pct, one probe per symbol
.sym <- P$GPL14550$sym; .md <- rowMeans(P$GPL14550$e[, main$sample_id[main$platform == "GPL14550"]])
.keep <- !is.na(.sym) & .sym != "" & .md > quantile(.md, 0.25); .ord <- order(-.md); .kk <- .ord[.keep[.ord]]; .kk <- .kk[!duplicated(.sym[.kk])]
genes <- rownames(P$GPL14550$e)[.kk]; gene_sym <- .sym[.kk]
FORM <- ~ DLCO + FEV1 + sqrtEmph + Age + Male + Current
fitD <- function(E, s, form = FORM) { X <- model.matrix(form, s); f <- eBayes(lmFit(E, X), trend = TRUE, robust = TRUE)
  list(t = f$t[, "DLCO"], p = f$p.value[, "DLCO"], b = f$coefficients[, "DLCO"], se = f$stdev.unscaled[, "DLCO"] * sqrt(f$s2.post), X = X) }
# Array QC (identical rule on both platforms): median inter-array correlation below Q1 - 3*IQR (Tukey far-out)
qc_outliers <- function(p) { s <- main[main$platform == p, ]; E <- P[[p]]$e[genes, s$sample_id]
  cc <- cor(E); mc <- apply(cc, 1, function(x) median(x[x < 1])); q <- quantile(mc, c(.25, .75))
  names(mc)[mc < q[1] - 3 * diff(q)] }
main$qc_outlier <- main$sample_id %in% c(qc_outliers("GPL14550"), qc_outliers("GPL6480"))
