# ============================================================
# Step 6: 决定性诊断 —— 多基因集富集 + 随机零分布 Z 分数
#  用已知阳性对照(星形胶质/小胶质/神经元标志物)检验流程有效性
#  若阳性对照强富集而乳酸化 Z≈0 -> 结论可靠: 乳酸化确无特异信号
# ============================================================
suppressMessages({library(limma); library(jsonlite)})
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); OUT <- file.path(BASE, "out")

expr_df <- read.delim(file.path(OUT, "expr_matrix.tsv"), row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
expr <- as.matrix(expr_df); storage.mode(expr) <- "double"
meta <- read.csv(file.path(OUT, "design_aligned.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
rownames(meta) <- meta$matrix_col; meta <- meta[colnames(expr), ]

map <- fromJSON(file.path(OUT, "ensg2symbol.json"))
sym <- sapply(sub("\\..*$", "", rownames(expr)), function(g) {
  v <- map[[g]]; if (is.null(v)) NA_character_ else if (is.list(v)) v$symbol[1] else NA_character_ })
names(sym) <- rownames(expr)
ensg <- sub("\\..*$", "", rownames(expr))

## ---------- limma ----------
keep <- meta$disease %in% c("TLE-HS", "ControlCortex")
e <- expr[, keep]; m <- meta[keep, ]
m$grp <- factor(ifelse(m$disease == "TLE-HS", "Case", "Ctrl"), levels = c("Ctrl", "Case"))
design <- model.matrix(~ grp + factor(tissue), data = m)
fit <- eBayes(lmFit(e, design), trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
tt$ensg <- ensg; tt$symbol <- sym[rownames(tt)]
is_deg <- tt$adj.P.Val < 0.05
cat(sprintf("[对比] TLE-HS(64) vs ControlCortex(14)   DEG=%.1f%%\n\n", 100*mean(is_deg)))

## ---------- 基因集 ----------
GS <- list(
 "★ 乳酸化 LRG (332)" = readLines(file.path(OUT, "LRG_final_symbols.txt")),
 "阳性对照: 星形胶质细胞" = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","GJA1","S100B",
                        "ALDOC","GLUL","SOX9","NFIA","NFIB","CD44","VIM"),
 "阳性对照: 小胶质细胞"   = c("TMEM119","P2RY12","CX3CR1","AIF1","CD68","CSF1R","TREM2",
                        "TYROBP","C1QA","C1QB","C1QC","ITGAM","SPI1","CD14"),
 "阳性对照: 神经元"       = c("RBFOX3","SYP","SNAP25","NEFL","NEFM","MAP2","TUBB3","GRIN1",
                        "SLC17A7","GAD1","GAD2","SYT1","NEFH","STMN2","UCHL1"),
 "阳性对照: 少突胶质细胞" = c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","SOX10","OLIG1","OLIG2"),
 "阳性对照: 内皮细胞"     = c("CLDN5","PECAM1","VWF","CDH5","FLT1","KDR","TEK","ENG","ESAM"),
 "其它热点: 铁死亡"       = c("GPX4","ACSL4","SLC7A11","FTH1","FTL","NCOA4","TFRC","ALOX15",
                        "SAT1","GCLC","GCLM","SLC3A2","PCBP1","PCBP2","HMOX1"),
 "其它热点: 铜死亡"       = c("FDX1","LIAS","LIPT1","DLD","DLAT","PDHA1","PDHB","MTF1","GLS","CDKN2A"),
 "通路: 糖酵解(KEGG)"     = c("HK1","HK2","HK3","GPI","PFKL","PFKM","PFKP","ALDOA","ALDOB","ALDOC",
                        "TPI1","GAPDH","PGK1","PGAM1","PGAM2","ENO1","ENO2","ENO3","PKM","PKLR",
                        "LDHA","LDHB","LDHC","PDHA1","PDHB","PC","ACLY"),
 "通路: 氧化磷酸化核心"   = c("NDUFA1","NDUFB1","NDUFS1","SDHA","SDHB","UQCRC1","UQCRC2","CYC1",
                        "COX4I1","COX5A","ATP5F1A","ATP5F1B","ATP5MC1","VDAC1","VDAC2"),
 "通路: 细胞凋亡核心"     = c("BAX","BAK1","CASP3","CASP7","CASP8","CASP9","BCL2","BCL2L1","MCL1",
                        "APAF1","CYCS","DIABLO","XIAP","BIRC5","FAS","FADD","TP53","PMAIP1")
)

## ---------- 随机零分布 ----------
set.seed(20260922); N <- 2000
null_cache <- list()
get_null <- function(k) {
  key <- as.character(k)
  if (!is.null(null_cache[[key]])) return(null_cache[[key]])
  v <- replicate(N, sum(is_deg[sample.int(nrow(tt), k)]))
  null_cache[[key]] <<- v; v
}

res <- data.frame()
for (nm in names(GS)) {
  g <- unique(GS[[nm]])
  hit <- ensg %in% { l <- unique(g); l }
  # 用 symbol 匹配（更稳）
  hit <- toupper(tt$symbol) %in% toupper(g)
  k <- sum(hit)
  if (k < 5) { cat(sprintf("%-26s 检出=%d  (跳过)\n", nm, k)); next }
  obs <- sum(is_deg[hit])
  nl <- get_null(k)
  z <- (obs - mean(nl)) / sd(nl)
  p <- (sum(nl >= obs) + 1) / (N + 1)
  frac_gs <- obs / k; frac_all <- mean(is_deg)
  res <- rbind(res, data.frame(gene_set=nm, n_tested=k, n_sig=obs,
                               pct_sig=round(100*frac_gs,1), null_mean=round(mean(nl),1),
                               Z=round(z,2), empP=signif(p,3)))
}
res <- res[order(-res$Z), ]
cat("===============================================================\n")
cat(" 多基因集富集分析 (按 Z 分数排序)\n")
cat("   Z>0 且经验 P<0.05 = 该基因集显著富集于差异表达基因\n")
cat("===============================================================\n")
print(res, row.names = FALSE)
cat(sprintf("\n 全基因组 DEG 基线 = %.1f%%\n", 100*mean(is_deg)))
write.csv(res, file.path(OUT, "geneset_enrichment_Z.csv"), row.names = FALSE)
cat("\n[保存] geneset_enrichment_Z.csv\n")
