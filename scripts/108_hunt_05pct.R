# Step 108: 穷举搜索能产生 ~0.5% DEG 率的模型
# 手稿称 "tissue-matched analysis ... fell from 44.3% to 0.5% (-99%)"
# 需确认 0.5% 究竟出自哪个模型，或根本不存在
suppressMessages({library(limma)})
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); OUT <- file.path(BASE, "out")

expr_df <- read.delim(file.path(OUT, "expr_matrix.tsv"), row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
expr <- as.matrix(expr_df); storage.mode(expr) <- "double"
meta <- read.csv(file.path(OUT, "design_aligned.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$matrix_col
meta <- meta[colnames(expr), ]

# 细胞得分（与 13_cell_adjusted.R 相同做法：标志基因 z 均值）
MK <- list(astro = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","GJA1","S100B","ALDOC","GLUL","SOX9","NFIA","NFIB","CD44","VIM"),
           micro = c("TMEM119","P2RY12","CX3CR1","AIF1","CD68","CSF1R","TREM2","TYROBP","C1QA","C1QB","C1QC","ITGAM","SPI1","CD14"),
           oligo = c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","SOX10","OLIG1","OLIG2"),
           neuro = c("RBFOX3","SYP","SNAP25","NEFL","NEFM","MAP2","TUBB3","GRIN1","SLC17A7","GAD1","GAD2","SYT1","NEFH","STMN2","UCHL1"),
           endo  = c("CLDN5","PECAM1","VWF","CDH5","FLT1","KDR","TEK","ENG","ESAM"),
           immune= c("PTPRC","CD3D","CD3E","CD2","CD8A","CD4","MS4A1","CD79A","CD79B","IGHM","IGKC","CD19","CD37","CD48","CORO1A"))
map <- jsonlite::fromJSON(file.path(OUT, "ensg2symbol.json"))
sym <- sapply(sub("\\..*$", "", rownames(expr)), function(g) {
  v <- map[[g]]; if (is.null(v)) NA_character_ else if (is.list(v)) v$symbol[1] else NA_character_ })
z <- t(scale(t(expr)))
scores <- sapply(MK, function(gs) {
  idx <- which(toupper(sym) %in% toupper(gs))
  if (!length(idx)) return(rep(NA_real_, ncol(expr)))
  colMeans(z[idx, , drop = FALSE], na.rm = TRUE)
})
scores <- as.data.frame(scores)
cat("细胞得分列:", paste(colnames(scores), collapse = ", "), "\n\n")

run <- function(keep, rhs, tag) {
  e <- expr[, keep]; m <- meta[keep, ]
  m$grp <- factor(ifelse(m$disease == "TLE-HS", "Case", "Ctrl"), levels = c("Ctrl", "Case"))
  m <- cbind(m, scores[keep, , drop = FALSE])
  f <- as.formula(paste("~ grp", rhs))
  design <- tryCatch(model.matrix(f, data = m), error = function(e) NULL)
  if (is.null(design)) { cat(sprintf("%-64s 模型失败\n", tag)); return(invisible()) }
  fit <- eBayes(lmFit(e, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
  n <- sum(tt$adj.P.Val < 0.05)
  cat(sprintf("%-64s n=%3d  DEG=%6d = %5.2f%%\n", tag, length(keep), n, 100 * n / nrow(tt)))
  invisible(100 * n / nrow(tt))
}

CC <- meta$disease %in% c("TLE-HS", "ControlCortex")
ALL <- meta$disease %in% c("TLE-HS", "ControlCortex", "ControlHippocampus")
TEMP <- meta$disease %in% c("TLE-HS", "ControlCortex") & meta$tissue == "Temporal"

cat("=== 候选模型 ===\n")
run(CC,   "",                                          "TLE-HS vs ControlCortex        无协变量")
run(CC,   " + factor(tissue)",                         "TLE-HS vs ControlCortex        + tissue")
run(CC,   " + astro+micro+oligo+neuro+endo+immune",    "TLE-HS vs ControlCortex        + 6 细胞得分")
run(CC,   " + factor(tissue) + astro+micro+oligo+neuro+endo+immune",
                                                       "TLE-HS vs ControlCortex        + tissue + 6 细胞得分")
run(ALL,  " + factor(tissue)",                         "TLE-HS vs 全部对照              + tissue")
run(ALL,  " + factor(tissue) + astro+micro+oligo+neuro+endo+immune",
                                                       "TLE-HS vs 全部对照              + tissue + 6 细胞得分")
run(ALL,  " + astro+micro+oligo+neuro+endo+immune",    "TLE-HS vs 全部对照              + 6 细胞得分")
cat("\n仅 Temporal 组织（Case 64 / Ctrl 4）:\n")
run(TEMP, "",                                          "  仅 Temporal                  无协变量")
run(TEMP, " + astro+micro+oligo+neuro+endo+immune",    "  仅 Temporal                  + 6 细胞得分")

cat("\n=== 全库扫描：所有 out/*.csv 中是否出现 0.5% 量级的 DEG 计数 ===\n")
fs <- list.files(OUT, pattern = "\\.csv$", full.names = TRUE)
for (f in fs) {
  d <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d)) next
  nm <- names(d)
  hit <- grep("deg|n_sig|n_deg", nm, ignore.case = TRUE, value = TRUE)
  if (!length(hit)) next
  for (h in hit) {
    v <- suppressWarnings(as.numeric(d[[h]]))
    if (any(!is.na(v) & v > 0 & v / 19999 < 0.01 & v / 19999 > 0.0001))
      cat(sprintf("  %-42s 列=%-14s 值=%s\n", basename(f), h,
                  paste(head(v[!is.na(v) & v > 0 & v / 19999 < 0.01], 5), collapse = ",")))
  }
}
