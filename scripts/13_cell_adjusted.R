# ============================================================
# Step 7: 细胞组成校正后重测 —— 决定项目生死的实验
#  假设: 44% DEG 由胶质/免疫细胞组成差异驱动, 淹没了真实的疾病信号
#  做法: 用 marker 基因构建细胞类型评分 -> 作为协变量 -> 重跑 limma -> 再测富集
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
names(sym) <- rownames(expr); ensg <- sub("\\..*$", "", rownames(expr))

MARK <- list(
 astro = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","GJA1","S100B","ALDOC","GLUL","SOX9","NFIA","NFIB","CD44","VIM"),
 micro = c("TMEM119","P2RY12","CX3CR1","AIF1","CD68","CSF1R","TREM2","TYROBP","C1QA","C1QB","C1QC","ITGAM","SPI1","CD14"),
 oligo = c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","SOX10","OLIG1","OLIG2"),
 neuro = c("RBFOX3","SYP","SNAP25","NEFL","NEFM","MAP2","TUBB3","GRIN1","SLC17A7","GAD1","GAD2","SYT1","NEFH","STMN2","UCHL1"),
 endo  = c("CLDN5","PECAM1","VWF","CDH5","FLT1","KDR","TEK","ENG","ESAM"),
 immune= c("CD3D","CD3E","CD2","CD8A","CD4","MS4A1","CD79A","CD79B","IGHM","CD37","SELPLG","IKZF1","WAS","MNDA","PYCARD")
)
# 每样本 z-score 后的平均 = 细胞评分
z <- t(scale(t(expr)))
for (nm in names(MARK)) {
  idx <- which(toupper(sym) %in% MARK[[nm]])
  cat(sprintf("  评分 %-7s 用 %d 个 marker\n", nm, length(idx)))
  meta[[paste0("sc_", nm)]] <- colMeans(z[idx, , drop = FALSE], na.rm = TRUE)
}

## ---------- 对比 ----------
keep <- meta$disease %in% c("TLE-HS", "ControlCortex")
e <- expr[, keep]; m <- meta[keep, ]
m$grp <- factor(ifelse(m$disease == "TLE-HS", "Case", "Ctrl"), levels = c("Ctrl", "Case"))

cat("\n=== 细胞评分在两组间的差异 (Wilcoxon) ===\n")
for (nm in names(MARK)) {
  v <- m[[paste0("sc_", nm)]]
  w <- wilcox.test(v ~ m$grp)
  cat(sprintf("  %-7s  Case=%.2f  Ctrl=%.2f   P=%.3g\n", nm,
              mean(v[m$grp=="Case"]), mean(v[m$grp=="Ctrl"]), w$p.value))
}

run <- function(form, tag) {
  design <- model.matrix(form, data = m)
  fit <- eBayes(lmFit(e, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
  tt$ensg <- ensg; tt$symbol <- sym[rownames(tt)]
  is_deg <- tt$adj.P.Val < 0.05
  cat(sprintf("\n[%s]\n  总基因=%d  DEG=%.1f%%\n", tag, nrow(tt), 100*mean(is_deg)))
  tt
}

## 模型 A: 原模型
ttA <- run(~ grp + factor(tissue), "模型A: ~grp + tissue (原)")

## 模型 B: 细胞组成校正
ttB <- run(~ grp + factor(tissue) + sc_astro + sc_micro + sc_oligo + sc_neuro + sc_endo + sc_immune,
           "模型B: ~grp + tissue + 6个细胞评分 (组成校正)")

## ---------- 富集 Z ----------
GS <- list(
 "★ 乳酸化 LRG"     = readLines(file.path(OUT, "LRG_final_symbols.txt")),
 "小胶质细胞 markers" = MARK$micro, "星形胶质 markers" = MARK$astro,
 "少突胶质 markers"   = MARK$oligo, "神经元 markers" = MARK$neuro,
 "铁死亡"            = c("GPX4","ACSL4","SLC7A11","FTH1","FTL","NCOA4","TFRC","ALOX15","SAT1","GCLC","GCLM","SLC3A2","PCBP1","PCBP2","HMOX1"),
 "糖酵解"            = c("HK1","HK2","HK3","GPI","PFKL","PFKM","PFKP","ALDOA","ALDOB","ALDOC","TPI1","GAPDH","PGK1","PGAM1","PGAM2","ENO1","ENO2","ENO3","PKM","PKLR","LDHA","LDHB","LDHC")
)
set.seed(20260922); N <- 2000
setZ <- function(tt, tag) {
  is_deg <- tt$adj.P.Val < 0.05
  cat(sprintf("\n=== %s ===\n", tag))
  out <- data.frame()
  for (nm in names(GS)) {
    hit <- toupper(tt$symbol) %in% toupper(unique(GS[[nm]]))
    k <- sum(hit); if (k < 5) next
    obs <- sum(is_deg[hit])
    nl <- replicate(N, sum(is_deg[sample.int(nrow(tt), k)]))
    zz <- (obs - mean(nl)) / sd(nl)
    p <- (sum(nl >= obs) + 1) / (N + 1)
    out <- rbind(out, data.frame(gene_set=nm, n=k, n_sig=obs,
                Z=round(zz,2), empP=signif(p,3)))
  }
  print(out[order(-out$Z), ], row.names = FALSE)
  out
}
rA <- setZ(ttA, "模型A (未校正)")
rB <- setZ(ttB, "模型B (细胞组成校正后)")

write.csv(rA, file.path(OUT,"enrich_modelA.csv"), row.names=FALSE)
write.csv(rB, file.path(OUT,"enrich_modelB_celladjusted.csv"), row.names=FALSE)
cat("\n[保存] enrich_modelA.csv / enrich_modelB_celladjusted.csv\n")
