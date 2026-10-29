# Step 117: 在"基因集分析口径"(有 HGNC symbol) 上统一重算 TableS4 与 Table 1 的相关数字
# 消除 11,905(19958) 与 11,923(19999) 的混用
suppressMessages({library(limma); library(edgeR)})
set.seed(20260922)
NPERM <- 1000
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); MOUT <- file.path(BASE, "out", "multi")
OUT <- file.path(BASE, "out")

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

rows <- list(); t1rows <- list()
for (ds in c("GSE256068", "GSE213488", "GSE128300", "GSE272142")) {
  E <- as.matrix(read.delim(file.path(MOUT, paste0(ds, "_expr.tsv")), row.names = 1,
                            check.names = FALSE))
  dd <- read.csv(file.path(MOUT, paste0(ds, "_design.csv")), stringsAsFactors = FALSE)
  dd <- dd[dd$sample %in% colnames(E), ]; E <- E[, dd$sample, drop = FALSE]
  E <- E[rowSums(is.finite(E)) == ncol(E), , drop = FALSE]

  ## 先统一到"有 symbol"口径，再做计数过滤/建 DGEList（否则 dge 与 E 行数不一致）
  sy <- E2S[sub("\\..*$", "", rownames(E))]
  k <- !is.na(sy); E <- E[k, , drop = FALSE]; sy <- sy[k]

  is_counts <- max(E, na.rm = TRUE) > 100
  if (is_counts) {
    keep <- rowSums(E >= 5) >= max(3, 0.2 * ncol(E)); E <- E[keep, , drop = FALSE]
    sy <- sy[keep]
    dge <- DGEList(counts = round(E)); dge <- calcNormFactors(dge)
  }
  grp <- factor(dd$group, levels = c("Control", "Case"))

  prep <- function(des) {
    if (is_counts) eBayes(lmFit(voom(dge, des, plot = FALSE), des), trend = TRUE, robust = TRUE)
    else eBayes(lmFit(E, des), trend = TRUE, robust = TRUE)
  }
  dA <- model.matrix(~ grp)
  ttA <- topTable(prep(dA), coef = "grpCase", number = Inf, sort.by = "P")
  ttA$symbol <- sy
  isA <- ttA$adj.P.Val < 0.05
  degA <- sum(isA, na.rm = TRUE)

  Z <- t(scale(t(E))); sc <- data.frame(row.names = colnames(E))
  for (nm in names(MARK)) {
    idx <- which(toupper(sy) %in% MARK[[nm]])
    if (length(idx) < 3) next
    sc[[nm]] <- colMeans(Z[idx, , drop = FALSE], na.rm = TRUE)
  }
  dB <- model.matrix(as.formula(paste("~ grp +", paste(names(sc), collapse = " + "))),
                     data = cbind(data.frame(grp = grp), sc))
  ttB <- topTable(prep(dB), coef = "grpCase", number = Inf, sort.by = "P")
  degB <- sum(ttB$adj.P.Val < 0.05, na.rm = TRUE)

  cat(sprintf("%-11s genes=%5d  ModelA=%6d (%.2f%%)  ModelB=%6d (%.2f%%)  降幅=%.1f%%\n",
              ds, nrow(ttA), degA, 100 * degA / nrow(ttA), degB, 100 * degB / nrow(ttB),
              100 * (degA - degB) / degA))

  for (nm in names(MARK)) {
    hit <- toupper(ttA$symbol) %in% MARK[[nm]]; kk <- sum(hit)
    if (kk < 5) { rows[[length(rows) + 1]] <- data.frame(dataset = ds, marker = nm, n = kk,
                     obs = NA, exp = NA, Z = NA); next }
    obs <- sum(isA[hit])
    br <- unique(quantile(ttA$AveExpr, probs = seq(0, 1, 0.05), na.rm = TRUE))
    st <- cut(ttA$AveExpr, breaks = br, include.lowest = TRUE, labels = FALSE); st[is.na(st)] <- 1L
    need <- tabulate(st[hit], nbins = max(st)); ib <- split(seq_along(st), st)
    sims <- replicate(NPERM, { t <- 0L
      for (s in seq_along(need)) { m <- need[s]; if (!m) next
        pool <- ib[[as.character(s)]]; if (is.null(pool) || !length(pool)) next
        t <- t + sum(isA[if (length(pool) <= m) pool else sample(pool, m)]) }; t })
    z <- if (sd(sims) > 0) (obs - mean(sims)) / sd(sims) else NA
    rows[[length(rows) + 1]] <- data.frame(dataset = ds, marker = nm, n = kk, obs = obs,
                                           exp = round(mean(sims), 1), Z = round(z, 2))
  }
  rows[[length(rows) + 1]] <- data.frame(dataset = ds, marker = "DEG_modelA",
                                         n = nrow(ttA), obs = degA, exp = NA, Z = NA)
  rows[[length(rows) + 1]] <- data.frame(dataset = ds, marker = "DEG_modelB_celladj",
                                         n = nrow(ttB), obs = degB, exp = NA, Z = NA)
  t1rows[[ds]] <- data.frame(dataset = ds, n_genes = nrow(ttA), degA = degA,
                             degB = degB, pctA = round(100 * degA / nrow(ttA), 2),
                             red = round(100 * (degA - degB) / degA, 1))
}
cmp <- do.call(rbind, rows)
write.csv(cmp, file.path(OUT, "cell_composition_analysis_symbolframe.csv"), row.names = FALSE)
t1 <- do.call(rbind, t1rows)
write.csv(t1, file.path(OUT, "table1_deg_symbolframe.csv"), row.names = FALSE)
cat("\n已写入 out/cell_composition_analysis_symbolframe.csv 与 out/table1_deg_symbolframe.csv\n")
print(t1, row.names = FALSE)
