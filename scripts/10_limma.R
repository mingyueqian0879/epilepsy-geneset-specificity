# ============================================================
# Step 4: limma 差异表达分析  —— TLE-HS vs Control
#  数据: GSE256068 (log2 TPM, RNA-seq)
#  对比: 主分析 TLE-HS vs ControlCortex ; 敏感性 TLE-HS vs All Control
#  核心问题: 332 个 LRG 中有多少显著差异表达?
# ============================================================
suppressMessages({
  library(limma)
  library(jsonlite)
})

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT  <- file.path(BASE, "out")

cat("========================================================\n")
cat(" limma:  TLE-HS vs Control   (GSE256068)\n")
cat("========================================================\n\n")

## ---------- 1. 读数据 ----------
cat("[读取表达矩阵] ...\n")
expr_df <- read.delim(file.path(OUT, "expr_matrix.tsv"), row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
expr <- as.matrix(expr_df)
storage.mode(expr) <- "double"
meta <- read.csv(file.path(OUT, "design_aligned.csv"), stringsAsFactors = FALSE,
                 fileEncoding = "UTF-8")
rownames(meta) <- meta$matrix_col
meta <- meta[colnames(expr), ]

cat("[数据] 基因数 =", nrow(expr), " 样本数 =", ncol(expr), "\n")
cat("[表达值] 范围", round(min(expr), 2), "~", round(max(expr), 2), "\n\n")

## ---------- 2. LRG ----------
lrg_ensg <- readLines(file.path(OUT, "LRG_final_ensg.txt"))
lrg_ensg <- unique(trimws(lrg_ensg))
lrg_ensg <- sub("\\..*$", "", lrg_ensg)
expr_ensg <- sub("\\..*$", "", rownames(expr))

lrg_in <- intersect(lrg_ensg, expr_ensg)
cat("[LRG] 列表共", length(lrg_ensg), "个 ENSG; 在表达矩阵中检出", length(lrg_in), "个\n\n")

## ---------- 3. 通用 limma 函数 ----------
run_limma <- function(expr, meta, grp_col, case, ctrl, covar = NULL, tag = "") {
  keep <- meta[[grp_col]] %in% c(case, ctrl)
  e <- expr[, keep, drop = FALSE]
  m <- meta[keep, , drop = FALSE]
  m$grp <- factor(ifelse(m[[grp_col]] == case, "Case", "Ctrl"),
                  levels = c("Ctrl", "Case"))
  cat("------------------------------------------------------\n")
  cat(sprintf("[%s]  Case(%s)=%d   Ctrl=%d\n", tag, case, sum(m$grp == "Case"), sum(m$grp == "Ctrl")))
  print(table(m$grp, m$tissue))

  if (is.null(covar)) {
    design <- model.matrix(~ grp, data = m)
  } else {
    m$.cov <- factor(m[[covar]])
    design <- model.matrix(~ grp + .cov, data = m)
  }
  fit <- lmFit(e, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
  tt$ensg <- sub("\\..*$", "", rownames(tt))
  tt$is_lrg <- tt$ensg %in% lrg_ensg
  attr(tt, "n_case") <- sum(m$grp == "Case")
  attr(tt, "n_ctrl") <- sum(m$grp == "Ctrl")
  tt
}

## ---------- 4. 主分析 ----------
tt_main <- run_limma(expr, meta, "disease", "TLE-HS", "ControlCortex",
                     covar = "tissue", tag = "主分析 TLE-HS vs ControlCortex")

## ---------- 5. 敏感性分析 ----------
meta$ctrl_all <- ifelse(grepl("^Control", meta$disease), "Control", meta$disease)
tt_sens <- run_limma(expr, meta, "ctrl_all", "TLE-HS", "Control",
                     covar = "tissue", tag = "敏感性 TLE-HS vs 全部对照")

## ---------- 6. 统计 LRG ----------
summarise_tt <- function(tt, label) {
  tt$ensg <- sub("\\..*$", "", rownames(tt))
  tt$is_lrg <- tt$ensg %in% lrg_ensg
  padj <- tt$adj.P.Val; pval <- tt$P.Value
  deg  <- sum(padj < 0.05, na.rm = TRUE)
  deg1 <- sum(padj < 0.05 & abs(tt$logFC) > 1, na.rm = TRUE)
  lrg_tested <- sum(tt$is_lrg)
  lrg_sig    <- sum(tt$is_lrg & padj < 0.05, na.rm = TRUE)
  lrg_sig01  <- sum(tt$is_lrg & padj < 0.01, na.rm = TRUE)
  lrg_p005   <- sum(tt$is_lrg & pval < 0.05, na.rm = TRUE)

  # Fisher 精确检验: LRG 是否富集于 DEG
  a <- sum(tt$is_lrg & padj < 0.05); b <- sum(tt$is_lrg & padj >= 0.05)
  c <- sum(!tt$is_lrg & padj < 0.05); d <- sum(!tt$is_lrg & padj >= 0.05)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")

  cat("\n======================================================\n")
  cat(" ", label, "\n")
  cat("======================================================\n")
  cat(sprintf("  总基因数            : %d\n", nrow(tt)))
  cat(sprintf("  DEG (adj.P<0.05)    : %d\n", deg))
  cat(sprintf("  DEG (adj.P<0.05 & |logFC|>1): %d\n", deg1))
  cat(sprintf("  ----- LRG -----\n"))
  cat(sprintf("  LRG 检出(可检验)    : %d\n", lrg_tested))
  cat(sprintf("  LRG 显著 adj.P<0.05 : %d   <== 关键指标\n", lrg_sig))
  cat(sprintf("  LRG 显著 adj.P<0.01 : %d\n", lrg_sig01))
  cat(sprintf("  LRG 名义 P<0.05     : %d\n", lrg_p005))
  cat(sprintf("  LRG 富集 Fisher OR  : %.2f  (P = %.3g)\n", ft$estimate, ft$p.value))
  cat(sprintf("  显著 LRG 占比       : %.1f%%  vs 全基因组 %.1f%%\n",
              100 * lrg_sig / max(lrg_tested, 1), 100 * deg / nrow(tt)))

  sig <- tt[tt$is_lrg & tt$adj.P.Val < 0.05, c("ensg", "logFC", "AveExpr", "P.Value", "adj.P.Val")]
  sig <- sig[order(sig$adj.P.Val), ]
  if (nrow(sig) > 0) {
    cat("\n  --- 显著 LRG 明细 (按 adj.P 排序) ---\n")
    print(head(sig, 40), row.names = FALSE, digits = 3)
  }
  invisible(list(n=nrow(tt), deg=deg, lrg_tested=lrg_tested, lrg_sig=lrg_sig,
                 lrg_sig01=lrg_sig01, lrg_p005=lrg_p005, fisher_or=as.numeric(ft$estimate),
                 fisher_p=ft$p.value, sig_table=sig))
}

res_main <- summarise_tt(tt_main, "主分析: TLE-HS vs ControlCortex")
res_sens <- summarise_tt(tt_sens, "敏感性: TLE-HS vs 全部对照")

## ---------- 7. 保存 ----------
write.csv(tt_main, file.path(OUT, "limma_TLEvsControlCortex.csv"), row.names = FALSE)
write.csv(tt_sens, file.path(OUT, "limma_TLEvsAllControl.csv"), row.names = FALSE)
saveRDS(list(main = res_main, sens = res_sens), file.path(OUT, "limma_results.rds"))

cat("\n[保存] limma_TLEvsControlCortex.csv / limma_TLEvsAllControl.csv\n")
