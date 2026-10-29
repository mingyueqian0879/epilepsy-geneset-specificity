# Step 115: 用脚本 36 的原始标记基因与预处理，重算 GSE256068 的四级模型阶梯
#   M1 ~ grp
#   M2 ~ grp + tissue
#   M3 ~ grp + 6 细胞得分
#   M4 ~ grp + tissue + 6 细胞得分
# 两个比较集：(i) TLE-HS vs 全部对照(27)  (ii) TLE-HS vs ControlCortex(14)
suppressMessages({library(limma); library(edgeR)})
set.seed(20260922)
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); MOUT <- file.path(BASE, "out", "multi")

mapdf <- read.csv(file.path(MOUT, "ensg2symbol_gencode.csv"), stringsAsFactors = FALSE)
E2S <- setNames(mapdf$symbol, mapdf$ensg)

MARK <- list(
 astro = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","GJA1","S100B","ALDOC","GLUL","SOX9",
           "NFIA","NFIB","CD44","VIM"),
 micro = c("TMEM119","P2RY12","CX3CR1","AIF1","CD68","CSF1R","TREM2","TYROBP","C1QA","C1QB",
           "C1QC","ITGAM","SPI1","CD14"),
 oligo = c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","SOX10","OLIG1","OLIG2"),
 neuro = c("RBFOX3","SYP","SNAP25","NEFL","NEFM","MAP2","TUBB3","GRIN1","SLC17A7","GAD1",
           "GAD2","SYT1","NEFH","STMN2","UCHL1"),
 endo  = c("CLDN5","PECAM1","VWF","CDH5","FLT1","KDR","TEK","ENG","ESAM"),
 immune= c("CD3D","CD3E","CD2","CD8A","CD4","MS4A1","CD79A","CD79B","IGHM","CD37","SELPLG",
           "IKZF1","WAS","MNDA","PYCARD"))

ds <- "GSE256068"
E <- as.matrix(read.delim(file.path(MOUT, paste0(ds, "_expr.tsv")), row.names = 1,
                          check.names = FALSE))
dd <- read.csv(file.path(MOUT, paste0(ds, "_design.csv")), encoding = "UTF-8",
               stringsAsFactors = FALSE)
dd <- dd[dd$sample %in% colnames(E), ]
E <- E[, dd$sample, drop = FALSE]
is_counts <- max(E, na.rm = TRUE) > 100
cat("is_counts =", is_counts, "  dim =", paste(dim(E), collapse = " x "), "\n")
if (is_counts) {
  keep <- rowSums(E >= 5) >= max(3, 0.2 * ncol(E)); E <- E[keep, , drop = FALSE]
  dge <- DGEList(counts = round(E)); dge <- calcNormFactors(dge)
} else {
  E <- E[rowSums(is.finite(E)) == ncol(E), , drop = FALSE]
}
cat("过滤后基因数 =", nrow(E), "\n\n")

Z <- t(scale(t(E)))
sc <- data.frame(row.names = colnames(E))
for (nm in names(MARK)) {
  idx <- which(toupper(E2S[sub("\\..*$", "", rownames(E))]) %in% MARK[[nm]])
  if (length(idx) < 3) next
  sc[[nm]] <- colMeans(Z[idx, , drop = FALSE], na.rm = TRUE)
}
cat("细胞得分列:", paste(names(sc), collapse = ", "), "\n\n")

run <- function(sel, form, tag) {
  ee <- E[, sel, drop = FALSE]; mm <- dd[sel, ]
  mm$grp <- factor(mm$group, levels = c("Control", "Case"))
  mm <- cbind(mm, sc[sel, , drop = FALSE])
  design <- model.matrix(as.formula(form), data = mm)
  if (is_counts) {
    vm <- voom(dge[, sel, drop = FALSE], design, plot = FALSE)
    fit <- eBayes(lmFit(vm, design), trend = TRUE, robust = TRUE)
  } else {
    fit <- eBayes(lmFit(ee, design), trend = TRUE, robust = TRUE)
  }
  tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
  n <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("  %-38s DEG=%6d / %5d = %6.2f%%\n", tag, n, nrow(tt), 100 * n / nrow(tt)))
  invisible(100 * n / nrow(tt))
}

ALL <- dd$detail %in% c("TLE-HS", "ControlCortex", "ControlHippocampus")
CC  <- dd$detail %in% c("TLE-HS", "ControlCortex")

cat("=== (i) TLE-HS(64) vs 全部对照(27) ===\n")
a1 <- run(ALL, "~ grp", "M1 仅分组")
a2 <- run(ALL, "~ grp + factor(tissue)", "M2 + 组织来源")
a3 <- run(ALL, "~ grp + astro+micro+oligo+neuro+endo+immune", "M3 + 6 细胞得分")
a4 <- run(ALL, "~ grp + factor(tissue) + astro+micro+oligo+neuro+endo+immune",
          "M4 + 组织来源 + 6 细胞得分")

cat("\n=== (ii) TLE-HS(64) vs ControlCortex(14) ===\n")
b1 <- run(CC, "~ grp", "M1 仅分组")
b2 <- run(CC, "~ grp + factor(tissue)", "M2 + 组织来源")
b3 <- run(CC, "~ grp + astro+micro+oligo+neuro+endo+immune", "M3 + 6 细胞得分")
b4 <- run(CC, "~ grp + factor(tissue) + astro+micro+oligo+neuro+endo+immune",
          "M4 + 组织来源 + 6 细胞得分")

cat("\n=== 与已存档 Table 1 对照 ===\n")
cat("  Table 1 报 Model A = 11905 (59.65%), Model B(细胞校正) = 4367 (-63.3%)\n")
cat(sprintf("  本次 (i) M1 = %.2f%%, M3 = %.2f%%\n", a1, a3))
