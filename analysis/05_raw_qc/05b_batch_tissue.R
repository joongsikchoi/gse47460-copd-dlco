# Step 5b: (a) tissue-composition profile of QC outliers, (b) array batch (scan date) vs DLCO, (c) H1/H2 with batch adjustment.
# Usage: Rscript 05b_batch_tissue.R <prepared dir> <01_cohort dir> <gmt> <raw_fe_headers csv> <out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; gmtf <- args[3]; fe <- args[4]; out <- args[5]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
suppressPackageStartupMessages(library(GSVA))
h <- read.csv(fe); main$scan_date <- h$scan_date[match(main$sample_id, h$sample_id)]
# (a) tissue-content markers
mk <- list(skeletal_muscle = c("MYH2","ACTA1","CKM","MB","TNNT3","MYL1","DES"), mesothelium = c("ITLN1","MSLN","UPK3B","CALB2","WT1"),
  blood_RBC = c("HBB","HBA1","HBA2","ALAS2","HBD"), airway_secretory = c("BPIFB1","SCGB1A1","MUC5B","MSMB","SCGB3A1"),
  cilia = c("FOXJ1","DNAH5","CAPS","TPPP3","RSPH1"), alveolar_AT2 = c("SFTPC","SFTPB","SFTPA1","NAPSA","LAMP3"),
  alveolar_cap = c("CA4","SOSTDC1","HPGD","FCN3","EDNRB"), adipose = c("ADIPOQ","PLIN1","LEP","CIDEA","FABP4"),
  lymphoid = c("CD19","MS4A1","CD79A","CXCL13","POU2AF1"), acute_stress = c("FOS","FOSB","EGR1","ATF3","IL6"))
tis <- do.call(rbind, lapply(names(P), function(p) { s <- main[main$platform == p, ]; e <- P[[p]]$e[, s$sample_id]; sym <- P[[p]]$sym
  sc <- sapply(mk, function(g) colMeans(t(scale(t(e[match(g[g %in% sym], sym), , drop = FALSE]))), na.rm = TRUE))
  data.frame(sample_id = s$sample_id, platform = p, qc_outlier = s$qc_outlier, round(sc, 2)) }))
write.csv(tis, file.path(out, "tissue_content_scores_cohort186.csv"), row.names = FALSE)
# (b) batch vs DLCO / residual DLCO
bt <- do.call(rbind, lapply(names(P), function(p) { s <- main[main$platform == p & !main$qc_outlier, ]
  s$resDLCO <- resid(lm(DLCO ~ FEV1 + sqrtEmph + Age + Male + Current, s))
  data.frame(platform = p, n = nrow(s), n_batches = length(unique(s$scan_date)),
    p_DLCO_by_batch = anova(lm(DLCO ~ factor(scan_date), s))[1, 5], p_residDLCO_by_batch = anova(lm(resDLCO ~ factor(scan_date), s))[1, 5],
    p_FEV1_by_batch = anova(lm(FEV1 ~ factor(scan_date), s))[1, 5]) }))
write.csv(bt, file.path(out, "batch_vs_clinical.csv"), row.names = FALSE); print(bt)
# (c) H1/H2 with batch (scan date) as fixed effect, QC cohort
gmt <- strsplit(readLines(gmtf), "\t"); gs <- setNames(lapply(gmt, `[`, -(1:2)), sub("\\|.*", "", sapply(gmt, `[`, 1))); h2sets <- grep("^H2_", names(gs), value = TRUE)
FB <- update(FORM, ~ . + factor(scan_date))
res <- list(); h2b <- data.frame()
for (p in names(P)) { s <- main[main$platform == p & !main$qc_outlier, ]; E <- P[[p]]$e[genes, s$sample_id]
  res[[p]] <- list(base = fitD(E, s), batch = fitD(E, s, FB))
  Eg <- E; rownames(Eg) <- gene_sym; sc <- t(scale(t(gsva(gsvaParam(Eg, gs[h2sets], minSize = 10), verbose = FALSE)))); s2 <- cbind(s, t(sc))
  for (g in h2sets) { f <- summary(lm(update(FB, as.formula(paste(g, "~ ."))), s2))$coef["DLCO", ]
    h2b <- rbind(h2b, data.frame(platform = p, set = g, beta_per10_batchadj = f[1] * 10, se = f[2] * 10, p = f[4])) } }
D <- res$GPL14550; R <- res$GPL6480; top <- order(D$batch$p)[1:100]
h1b <- data.frame(analysis = c("QC, no batch", "QC + scan-date batch"),
  disc_FDR05 = c(sum(p.adjust(D$base$p, "BH") < .05), sum(p.adjust(D$batch$p, "BH") < .05)),
  rep_FDR05 = c(sum(p.adjust(R$base$p, "BH") < .05), sum(p.adjust(R$batch$p, "BH") < .05)),
  conc_top100 = c(mean(sign(D$base$t[order(D$base$p)[1:100]]) == sign(R$base$t[order(D$base$p)[1:100]])), mean(sign(D$batch$t[top]) == sign(R$batch$t[top]))),
  tcor_disc_vs_nobatch = c(1, cor(D$base$t, D$batch$t)))
h2b$meta <- NA
for (g in h2sets) { d <- h2b[h2b$set == g, ]; w <- 1 / d$se^2; b <- sum(w * d$beta_per10_batchadj) / sum(w); se <- sqrt(1 / sum(w)); h2b$meta[h2b$set == g] <- sprintf("%.3f (p=%.4f)", b, 2 * pnorm(-abs(b / se))) }
write.csv(h1b, file.path(out, "H1_batch_adjusted.csv"), row.names = FALSE); write.csv(h2b, file.path(out, "H2_batch_adjusted.csv"), row.names = FALSE)
print(h1b, digits = 3); print(h2b, digits = 3)
