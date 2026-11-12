# Step 9b: camera repeated with inter-gene correlation ESTIMATED from the data (conservative), replaces the fixed 0.01 default.
# Usage: Rscript 09b_camera_estimated_cor.R <prepared dir> <01_cohort dir> <gmt> <09 out dir>
args <- commandArgs(TRUE); prep <- args[1]; cdir <- args[2]; gmtf <- args[3]; out <- args[4]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(GO.db) })
SC <- readRDS(file.path(out, "SC_with_tertile.rds"))
gmt <- strsplit(readLines(gmtf), "\t"); gs <- setNames(lapply(gmt, `[`, -(1:2)), sub("\\|.*", "", sapply(gmt, `[`, 1)))
go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(gene_sym), keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL")))
go <- unique(go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", c("SYMBOL", "GOALL")])
idx <- lapply(split(match(go$SYMBOL, gene_sym), go$GOALL), function(i) unique(i[!is.na(i)])); idx <- idx[lengths(idx) >= 15 & lengths(idx) <= 300]
ct <- lapply(gs, function(g) which(gene_sym %in% g))
run <- function(sets) lapply(names(P), function(p) { s <- SC[[p]]
  camera(P[[p]]$e[genes, s$sample_id], sets, model.matrix(FORM, s), contrast = "DLCO", inter.gene.cor = NA, sort = FALSE) })
cg <- run(idx); cc <- run(ct)
F6 <- data.frame(GO = names(idx), term = suppressMessages(AnnotationDbi::select(GO.db, names(idx), "TERM", "GOID"))$TERM, size = lengths(idx),
  disc_dir = cg[[1]]$Direction, disc_p = cg[[1]]$PValue, disc_FDR = cg[[1]]$FDR, rep_dir = cg[[2]]$Direction, rep_p = cg[[2]]$PValue, rep_FDR = cg[[2]]$FDR)
F6$replicated <- F6$disc_FDR < .05 & F6$rep_p < .05 & F6$disc_dir == F6$rep_dir; F6 <- F6[order(F6$disc_p), ]
write.csv(F6, file.path(out, "F_camera_GO_BP.csv"), row.names = FALSE)
cat("GO BP:", nrow(F6), "disc FDR<.05:", sum(F6$disc_FDR < .05), "rep FDR<.05:", sum(F6$rep_FDR < .05), "replicated:", sum(F6$replicated),
    " min disc FDR:", signif(min(F6$disc_FDR), 3), "\n")
top <- head(F6, 200); cat("direction agreement, top-200 discovery GO terms:", mean(top$disc_dir == top$rep_dir), " rep p<.05 among them:", sum(top$rep_p < .05 & top$disc_dir == top$rep_dir), "\n")
print(head(F6[, c("term", "size", "disc_dir", "disc_p", "disc_FDR", "rep_dir", "rep_p")], 15), digits = 3)
ctab <- do.call(rbind, lapply(1:2, function(k) data.frame(platform = names(P)[k], set = names(ct), NGenes = cc[[k]]$NGenes, Direction = cc[[k]]$Direction,
  PValue = cc[[k]]$PValue, FDR = p.adjust(cc[[k]]$PValue, "BH"))))
write.csv(ctab, file.path(out, "F_camera_celltype_sets.csv"), row.names = FALSE); print(ctab, digits = 3)
