# Step 116: 在"基因集分析口径"(有 HGNC symbol, 19958 基因) 上重算完整模型阶梯
# 目的：让 Table 1、Fig 2A、正文 §3.2/§3.9 使用同一口径，消除 11905/11923 混用
suppressMessages({library(limma)})
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

E <- as.matrix(read.delim(file.path(MOUT, "GSE256068_expr.tsv"), row.names = 1,
                          check.names = FALSE))
dd <- read.csv(file.path(MOUT, "GSE256068_design.csv"), stringsAsFactors = FALSE)
dd <- dd[dd$sample %in% colnames(E), ]; E <- E[, dd$sample, drop = FALSE]
E <- E[rowSums(is.finite(E)) == ncol(E), , drop = FALSE]

sym <- E2S[sub("\\..*$", "", rownames(E))]
cat("全矩阵基因数        :", nrow(E), " (有 symbol:", sum(!is.na(sym)), ")\n")

run_all <- function(sel, tag) {
  out <- list()
  for (frame in c("full", "symbol")) {
    if (frame == "full") { Ef <- E; sf <- sym } else {
      k <- !is.na(sym); Ef <- E[k, , drop = FALSE]; sf <- sym[k]
    }
    Z <- t(scale(t(Ef)))
    sc <- data.frame(row.names = colnames(Ef))
    for (nm in names(MARK)) {
      idx <- which(toupper(sf) %in% MARK[[nm]])
      if (length(idx) < 3) next
      sc[[nm]] <- colMeans(Z[idx, , drop = FALSE], na.rm = TRUE)
    }
    for (m in c("M1", "M2", "M3", "M4")) {
      ee <- Ef[, sel, drop = FALSE]; mm <- dd[sel, ]
      mm$grp <- factor(mm$group, levels = c("Control", "Case"))
      mm <- cbind(mm, sc[sel, , drop = FALSE])
      f <- switch(m, M1 = "~ grp", M2 = "~ grp + factor(tissue)",
                  M3 = "~ grp + astro+micro+oligo+neuro+endo+immune",
                  M4 = "~ grp + factor(tissue) + astro+micro+oligo+neuro+endo+immune")
      des <- model.matrix(as.formula(f), data = mm)
      tt <- topTable(eBayes(lmFit(ee, des), trend = TRUE, robust = TRUE),
                     coef = "grpCase", number = Inf, sort.by = "P")
      n <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
      out[[paste(frame, m)]] <- c(n = n, tot = nrow(tt), pct = 100 * n / nrow(tt))
    }
  }
  cat("\n---", tag, "---\n")
  for (frame in c("full", "symbol")) {
    lab <- if (frame == "full") "全矩阵(19999)" else "基因集口径(有symbol)"
    cat(sprintf("  %-22s", lab))
    for (m in c("M1", "M2", "M3", "M4")) {
      v <- out[[paste(frame, m)]]
      cat(sprintf("  %s=%5d(%5.2f%%)", m, v["n"], v["pct"]))
    }
    cat("\n")
  }
  invisible(out)
}

ALL <- dd$detail %in% c("TLE-HS", "ControlCortex", "ControlHippocampus")
CC  <- dd$detail %in% c("TLE-HS", "ControlCortex")
a <- run_all(ALL, "(i) TLE-HS(64) vs 全部对照(27)")
b <- run_all(CC,  "(ii) TLE-HS(64) vs ControlCortex(14)")

cat("\n=== 建议写入手稿的口径（基因集分析集，有 symbol）===\n")
cat(sprintf("  (i)  M1=%.1f%%  M2=%.1f%%  M3=%.1f%%  M4=%.2f%%\n",
            a[["symbol M1"]]["pct"], a[["symbol M2"]]["pct"],
            a[["symbol M3"]]["pct"], a[["symbol M4"]]["pct"]))
cat(sprintf("  (ii) M1=%.1f%%  M2=%.1f%%  M3=%.1f%%  M4=%.2f%%\n",
            b[["symbol M1"]]["pct"], b[["symbol M2"]]["pct"],
            b[["symbol M3"]]["pct"], b[["symbol M4"]]["pct"]))
json <- sprintf('{"all":{"M1":%.2f,"M2":%.2f,"M3":%.2f,"M4":%.2f,"M1n":%d,"M3n":%d},'
                '"cortex":{"M1":%.2f,"M2":%.2f,"M3":%.2f,"M4":%.2f,"M1n":%d,"M3n":%d}}',
                a[["symbol M1"]]["pct"], a[["symbol M2"]]["pct"], a[["symbol M3"]]["pct"],
                a[["symbol M4"]]["pct"], a[["symbol M1"]]["n"], a[["symbol M3"]]["n"],
                b[["symbol M1"]]["pct"], b[["symbol M2"]]["pct"], b[["symbol M3"]]["pct"],
                b[["symbol M4"]]["pct"], b[["symbol M1"]]["n"], b[["symbol M3"]]["n"])
writeLines(json, file.path(BASE, "out", "P0_model_ladder.json"))
cat("\n已写入 out/P0_model_ladder.json\n")
