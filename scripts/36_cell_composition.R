# ============================================================
# Step 36: 细胞组成混杂的多数据集量化
#  对每个数据集: DEG(原模型) vs DEG(+细胞评分校正)
#  以及: 细胞类型标志物的 Z 分数 (证明"信号=细胞组成")
# ============================================================
suppressMessages({library(limma); library(edgeR)})
set.seed(20260922)
NPERM <- 1000

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
MOUT <- file.path(BASE, "out", "multi")
OUT  <- file.path(BASE, "out")

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

DATASETS <- c("GSE256068", "GSE213488", "GSE128300", "GSE272142")
rows <- list()

for (ds in DATASETS) {
  E <- as.matrix(read.delim(file.path(MOUT, paste0(ds, "_expr.tsv")),
                            row.names = 1, check.names = FALSE))
  dd <- read.csv(file.path(MOUT, paste0(ds, "_design.csv")), encoding = "UTF-8",
                 stringsAsFactors = FALSE)
  dd <- dd[dd$sample %in% colnames(E), ]
  E <- E[, dd$sample, drop = FALSE]
  is_counts <- max(E, na.rm = TRUE) > 100
  if (is_counts) {
    keep <- rowSums(E >= 5) >= max(3, 0.2 * ncol(E))
    E <- E[keep, , drop = FALSE]
    dge <- DGEList(counts = round(E)); dge <- calcNormFactors(dge)
  } else {
    E <- E[rowSums(is.finite(E)) == ncol(E), , drop = FALSE]
  }
  grp <- factor(dd$group, levels = c("Control", "Case"))
  cat("\n", strrep("=", 74), "\n [", ds, "]  n=", ncol(E), "\n", sep = "")

  prep <- function(design) {
    if (is_counts) {
      vm <- voom(dge, design, plot = FALSE)
      eBayes(lmFit(vm, design), trend = TRUE, robust = TRUE)
    } else eBayes(lmFit(E, design), trend = TRUE, robust = TRUE)
  }

  ## 模型 A: 仅分组
  dA <- model.matrix(~ grp)
  ttA <- topTable(prep(dA), coef = "grpCase", number = Inf, sort.by = "P")
  ttA$symbol <- E2S[sub("\\..*$", "", rownames(ttA))]
  degA <- sum(ttA$adj.P.Val < 0.05, na.rm = TRUE)

  ## 细胞评分 (z-score 后取 marker 均值)
  Z <- t(scale(t(E)))
  sc <- data.frame(row.names = colnames(E))
  for (nm in names(MARK)) {
    idx <- which(toupper(E2S[sub("\\..*$", "", rownames(E))]) %in% MARK[[nm]])
    if (length(idx) < 3) next
    sc[[nm]] <- colMeans(Z[idx, , drop = FALSE], na.rm = TRUE)
  }
  ## 模型 B: 校正细胞评分
  formB <- as.formula(paste("~ grp +", paste(names(sc), collapse = " + ")))
  dB <- model.matrix(formB, data = cbind(data.frame(grp = grp), sc))
  ttB <- topTable(prep(dB), coef = "grpCase", number = Inf, sort.by = "P")
  degB <- sum(ttB$adj.P.Val < 0.05, na.rm = TRUE)

  cat(sprintf("  DEG: 模型A(仅分组) = %d (%.1f%%)   模型B(+%d细胞评分) = %d (%.1f%%)\n",
              degA, 100 * degA / nrow(ttA), ncol(sc), degB, 100 * degB / nrow(ttB)))

  ## 细胞评分在两组间的差异
  cat("  细胞评分组间差异 (Wilcoxon):\n")
  for (nm in colnames(sc)) {
    v <- sc[[nm]]
    p <- tryCatch(wilcox.test(v ~ grp)$p.value, error = function(e) NA)
    cat(sprintf("     %-7s Case=%+.2f Ctrl=%+.2f  P=%.2e\n", nm,
                mean(v[grp == "Case"]), mean(v[grp == "Control"]), p))
  }
  ## 阳性对照 Z (模型A下)
  isA <- ttA$adj.P.Val < 0.05
  cat("  细胞标志物基因集 Z (模型A):\n")
  for (nm in names(MARK)) {
    hit <- toupper(ttA$symbol) %in% MARK[[nm]]
    k <- sum(hit)
    if (k < 5) next
    obs <- sum(isA[hit])
    br <- unique(quantile(ttA$AveExpr, probs = seq(0, 1, 0.05), na.rm = TRUE))
    st <- cut(ttA$AveExpr, breaks = br, include.lowest = TRUE, labels = FALSE); st[is.na(st)] <- 1L
    need <- tabulate(st[hit], nbins = max(st)); ib <- split(seq_along(st), st)
    sims <- replicate(NPERM, {
      t <- 0L
      for (s in seq_along(need)) {
        m <- need[s]; if (!m) next
        pool <- ib[[as.character(s)]]; if (is.null(pool) || !length(pool)) next
        t <- t + sum(isA[if (length(pool) <= m) pool else sample(pool, m)])
      }; t })
    z <- if (sd(sims) > 0) (obs - mean(sims)) / sd(sims) else NA
    cat(sprintf("     %-7s n=%-3d obs=%-4d exp=%-7.1f Z=%+.2f\n", nm, k, obs, mean(sims), z))
    rows[[length(rows) + 1]] <- data.frame(dataset = ds, marker = nm, n = k, obs = obs,
                                           exp = round(mean(sims), 1), Z = round(z, 2))
  }
  rows[[length(rows) + 1]] <- data.frame(dataset = ds, marker = "DEG_modelA", n = nrow(ttA),
                                          obs = degA, exp = NA, Z = NA)
  rows[[length(rows) + 1]] <- data.frame(dataset = ds, marker = "DEG_modelB_celladj",
                                          n = nrow(ttB), obs = degB, exp = NA, Z = NA)
}

R <- do.call(rbind, rows)
write.csv(R, file.path(OUT, "cell_composition_analysis.csv"), row.names = FALSE)
cat("\n[saved] cell_composition_analysis.csv\n")
