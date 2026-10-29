# ============================================================
# Step 8: 挽救性验证 —— TLE-HS 内部: 乳酸化基因 vs 癫痫严重程度
#  设计优势: 全部为 TLE-HS 手术标本, 无"病例vs尸检"混杂
#            (细胞组成仍有变异, 但已用细胞评分校正)
#  GSE256068 提供 seisure frequency 字段
# ============================================================
suppressMessages({library(limma); library(jsonlite)})
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); OUT <- file.path(BASE, "out")

## 重新解析含 seizure frequency 的元数据
tags <- list()
con <- gzfile(file.path(BASE,"data","GSE256068_series_matrix.txt.gz"), "rt")
ln <- readLines(con, warn = FALSE); close(con)
for (l in ln) {
  if (startsWith(l, "!series_matrix_table_begin")) break
  if (!startsWith(l, "!Sample_")) next
  p <- strsplit(l, "\t")[[1]]
  tags[[p[1]]] <- rbind(tags[[p[1]]], sapply(p[-1], function(x) gsub('"','',x)))
}
n <- length(tags[["!Sample_geo_accession"]][1,])
grab <- function(key) {
  for (k in names(tags)) {
    if (!grepl("characteristics", k)) next
    m <- tags[[k]]
    for (i in seq_len(nrow(m))) {
      if (startsWith(tolower(m[i,1]), tolower(key)))
        return(sapply(m[i,], function(x) trimws(sub("^[^:]*:", "", x))))
    }
  }
  NULL
}
meta <- data.frame(
  matrix_col = as.character(tags[["!Sample_description"]][1,]),
  disease    = grab("disease state"),
  tissue     = grab("tissue"),
  age        = grab("age"),
  sex        = grab("Sex"),
  sez        = grab("seisure frequency"),
  stringsAsFactors = FALSE)
meta$sez_num <- suppressWarnings(as.numeric(meta$sez))
cat("[seizure frequency] 非缺失 =", sum(!is.na(meta$sez_num)),
    " / ", nrow(meta), "\n")
print(table(meta$disease, is.na(meta$sez_num)))

expr_df <- read.delim(file.path(OUT, "expr_matrix.tsv"), row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
expr <- as.matrix(expr_df); storage.mode(expr) <- "double"
rownames(meta) <- meta$matrix_col
meta <- meta[colnames(expr), ]
map <- fromJSON(file.path(OUT, "ensg2symbol.json"))
sym <- sapply(sub("\\..*$", "", rownames(expr)), function(g) {
  v <- map[[g]]; if (is.null(v)) NA_character_ else if (is.list(v)) v$symbol[1] else NA_character_ })
ensg <- sub("\\..*$", "", rownames(expr))

## ---------- TLE-HS 内部 ----------
sub <- meta$disease == "TLE-HS" & !is.na(meta$sez_num)
cat(sprintf("\n[TLE-HS 内部] n = %d (有发作频率)\n", sum(sub)))
e <- expr[, sub]; m <- meta[sub, ]
cat("  发作频率分布: "); print(summary(m$sez_num))

## 细胞评分(校正用)
z <- t(scale(t(e)))
MARK <- list(
 astro=c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","GJA1","S100B","ALDOC","GLUL","SOX9","NFIA","NFIB","CD44","VIM"),
 micro=c("TMEM119","P2RY12","CX3CR1","AIF1","CD68","CSF1R","TREM2","TYROBP","C1QA","C1QB","C1QC","ITGAM","SPI1","CD14"),
 oligo=c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","SOX10","OLIG1","OLIG2"),
 neuro=c("RBFOX3","SYP","SNAP25","NEFL","NEFM","MAP2","TUBB3","GRIN1","SLC17A7","GAD1","GAD2","SYT1","NEFH","STMN2","UCHL1"),
 immune=c("CD3D","CD3E","CD2","CD8A","CD4","MS4A1","CD79A","CD79B","IGHM","CD37","SELPLG","IKZF1","WAS","MNDA","PYCARD"))
for (nm in names(MARK)) {
  idx <- which(toupper(sym) %in% MARK[[nm]])
  m[[paste0("sc_",nm)]] <- colMeans(z[idx,,drop=FALSE], na.rm=TRUE)
}

## limma: 表达 ~ log2(发作频率) + 细胞评分   (连续变量关联)
m$sez_log <- log2(m$sez_num + 1)
design <- model.matrix(~ sez_log + sc_astro + sc_micro + sc_oligo + sc_neuro + sc_immune, data = m)
fit <- eBayes(lmFit(e, design), trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "sez_log", number = Inf, sort.by = "P")
tt$ensg <- ensg; tt$symbol <- sym[rownames(tt)]
is_sig <- tt$adj.P.Val < 0.05
cat(sprintf("\n[TLE内部 关联检验] 总基因=%d  显著(adj.P<0.05)=%d (%.2f%%)\n",
            nrow(tt), sum(is_sig), 100*mean(is_sig)))

## ---------- 富集 Z ----------
GS <- list(
 "★ 乳酸化 LRG"    = readLines(file.path(OUT,"LRG_final_symbols.txt")),
 "小胶质 markers"   = MARK$micro, "星形胶质 markers" = MARK$astro,
 "少突胶质 markers" = MARK$oligo, "神经元 markers" = MARK$neuro, "免疫 markers" = MARK$immune,
 "铁死亡" = c("GPX4","ACSL4","SLC7A11","FTH1","FTL","NCOA4","TFRC","ALOX15","SAT1","GCLC","GCLM","SLC3A2","PCBP1","PCBP2","HMOX1"),
 "糖酵解" = c("HK1","HK2","HK3","GPI","PFKL","PFKM","PFKP","ALDOA","ALDOB","ALDOC","TPI1","GAPDH","PGK1","PGAM1","PGAM2","ENO1","ENO2","ENO3","PKM","PKLR","LDHA","LDHB","LDHC"),
 "细胞凋亡"= c("BAX","BAK1","CASP3","CASP7","CASP8","CASP9","BCL2","BCL2L1","MCL1","APAF1","CYCS","DIABLO","XIAP","BIRC5","FAS","FADD","TP53","PMAIP1"))
set.seed(20260922); N <- 2000
res <- data.frame()
for (nm in names(GS)) {
  hit <- toupper(tt$symbol) %in% toupper(unique(GS[[nm]]))
  k <- sum(hit); if (k < 5) next
  obs <- sum(is_sig[hit])
  nl <- replicate(N, sum(is_sig[sample.int(nrow(tt), k)]))
  res <- rbind(res, data.frame(gene_set=nm, n=k, n_sig=obs,
     pct=round(100*obs/k,1), null=round(mean(nl),1),
     Z=round((obs-mean(nl))/sd(nl),2), empP=signif((sum(nl>=obs)+1)/(N+1),3)))
}
cat("\n===============================================================\n")
cat(" TLE-HS 内部: 乳酸化基因 vs 发作频率 (已校正细胞组成)\n")
cat("===============================================================\n")
print(res[order(-res$Z),], row.names=FALSE)
cat(sprintf("\n 全基因组显著基线 = %.2f%%\n", 100*mean(is_sig)))
write.csv(res, file.path(OUT,"enrich_TLE_internal_seizure.csv"), row.names=FALSE)
write.csv(tt,  file.path(OUT,"limma_TLE_internal_seizure.csv"), row.names=FALSE)

## 显著 LRG 明细
sl <- tt[toupper(tt$symbol) %in% toupper(readLines(file.path(OUT,"LRG_final_symbols.txt"))) & is_sig, ]
cat(sprintf("\n[TLE内部 显著 LRG] %d 个\n", nrow(sl)))
if (nrow(sl)>0) print(head(sl[,c("symbol","logFC","P.Value","adj.P.Val")],25), row.names=FALSE, digits=3)
cat("\n[保存] enrich_TLE_internal_seizure.csv / limma_TLE_internal_seizure.csv\n")
