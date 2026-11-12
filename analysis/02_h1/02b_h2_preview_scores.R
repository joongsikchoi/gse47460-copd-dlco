# EXPLORATORY preview (not the pre-specified H2 test): canonical alveolar marker scores vs DLCO.
# Gene sets are textbook markers chosen AFTER seeing discovery PC2 -> must not be reported as confirmatory H2.
# Usage: Rscript 02b_h2_preview_scores.R <prepared dir> <01_cohort dir> <out dir>
args <- commandArgs(TRUE); prep <- args[1]; d <- read.delim(file.path(args[2], "cohort_main_186.tsv")); out <- args[3]
d$cur <- d$Smoking == "Current"
sets <- list(aCap = c("CA4","SOSTDC1","HPGD","TBX2","EDNRB","APLN"), gCap = c("FCN3","GPIHBP1","IL7R","APLNR","PLVAP","VIPR1"),
  AT1 = c("AGER","RTKN2","CLIC5","EMP2","HOPX","GPRC5A"), AT2 = c("SFTPC","SFTPB","SFTPA1","SFTPA2","LAMP3","ABCA3","NAPSA","PGC"),
  stressIEG = c("FOSB","ATF3","IL6","EGR1","SOCS3","ZFP36"))
rows <- list()
for (p in c("GPL14550","GPL6480")) { s <- d[d$platform == p, ]
  e <- as.matrix(read.delim(gzfile(file.path(prep, paste0(p, "_expression.tsv.gz"))), row.names = 1, check.names = FALSE))[, s$sample_id]
  a <- read.delim(gzfile(file.path(prep, paste0(p, "_gene_annotation.tsv.gz"))), quote = "", colClasses = "character")
  sym <- a$GENE_SYMBOL[match(rownames(e), a$ID)]
  for (n in names(sets)) { g <- sets[[n]][sets[[n]] %in% sym]; sc <- colMeans(t(scale(t(e[match(g, sym), , drop = FALSE]))))
    f <- summary(lm(sc ~ DLCO + FEV1 + sqrt(Emphysema) + Age + Sex + cur, s))$coef["DLCO", ]
    rows[[length(rows) + 1]] <- data.frame(platform = p, set = n, n_genes = length(g),
      rho_DLCO_raw = cor(sc, s$DLCO, method = "s"), rho_FEV1_raw = cor(sc, s$FEV1, method = "s"), rho_emph_raw = cor(sc, s$Emphysema, method = "s"),
      adj_beta_per10 = f[1] * 10, adj_p = f[4]) } }
r <- do.call(rbind, rows); write.csv(r, file.path(out, "h2_preview_marker_scores.csv"), row.names = FALSE); print(r, digits = 3)
