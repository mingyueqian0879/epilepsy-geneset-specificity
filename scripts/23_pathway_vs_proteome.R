# ============================================================
# Step ②-E: 决定性对照 —— "蛋白质组衍生列表" vs "真实通路基因集"
#  假设: 乳酸化 LRG(332) 是【底物蛋白】列表, 功能高度异质 -> 无富集是必然
#        而离子通道/突触基因是【功能通路】, 应有信号
#  在两个独立层面检验: (a) bulk 转录组 Z 分数  (b) MR 命中率
# ============================================================
suppressMessages({library(limma); library(jsonlite)})
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); OUT <- file.path(BASE,"out")

## ---------- (a) bulk 转录组 ----------
expr_df <- read.delim(file.path(OUT,"expr_matrix.tsv"), row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
expr <- as.matrix(expr_df); storage.mode(expr) <- "double"
meta <- read.csv(file.path(OUT,"design_aligned.csv"), stringsAsFactors=FALSE, fileEncoding="UTF-8")
rownames(meta) <- meta$matrix_col; meta <- meta[colnames(expr),]
map <- fromJSON(file.path(OUT,"ensg2symbol.json"))
sym <- sapply(sub("\\..*$","",rownames(expr)), function(g){ v<-map[[g]]
  if(is.null(v)) NA_character_ else if(is.list(v)) v$symbol[1] else NA_character_})
ensg <- sub("\\..*$","",rownames(expr))

keep <- meta$disease %in% c("TLE-HS","ControlCortex")
e <- expr[,keep]; m <- meta[keep,]
m$grp <- factor(ifelse(m$disease=="TLE-HS","Case","Ctrl"), levels=c("Ctrl","Case"))
design <- model.matrix(~ grp + factor(tissue), data=m)
fit <- eBayes(lmFit(e,design), trend=TRUE, robust=TRUE)
tt <- topTable(fit, coef="grpCase", number=Inf, sort.by="P")
tt$symbol <- sym[rownames(tt)]; is_deg <- tt$adj.P.Val < 0.05
cat(sprintf("[bulk] TLE-HS(64) vs ControlCortex(14)  DEG=%.1f%%\n", 100*mean(is_deg)))

ION <- c("SCN1A","SCN1B","SCN2A","SCN3A","SCN8A","SCN9A","SCN1G","KCNQ2","KCNQ3","KCNT1","KCNT2",
 "KCNB1","KCNH1","KCNH5","KCNMA1","KCNC1","KCNA1","KCNA2","CACNA1A","CACNA1H","CACNA1E","CACNB4",
 "CACNA2D2","GABRA1","GABRB2","GABRB3","GABRG2","GABRD","GABBR2","GRIN1","GRIN2A","GRIN2B","GRIA3",
 "GRIA4","GRIK1","GRIK2","GRIK4","GRIK5","HCN1","HCN2","HCN4","CLCN2","CLCN4","ATP1A2","ATP1A3","ATP1B1")
SYN <- c("STXBP1","STX1A","STX1B","SYT1","SYT2","SYT11","SV2A","DNM1","DNM2","AP3B2","VAMP2",
 "SNAP25","PRRT2","UNC13A","UNC13B","RIMS1","RIMS2","NBEA","MUNC18","NSF","SNAP29","VAMP1")
LACREG <- c("EP300","CREBBP","KAT2A","KAT2B","KAT5","KAT7","KAT8","HDAC1","HDAC2","HDAC3",
 "HDAC8","HDAC10","HDAC11","SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7","ELP3")
GLYC <- c("HK1","HK2","HK3","GPI","PFKL","PFKM","PFKP","ALDOA","ALDOB","ALDOC","TPI1","GAPDH",
 "PGK1","PGAM1","PGAM2","ENO1","ENO2","ENO3","PKM","PKLR","LDHA","LDHB","LDHC")
LRG <- readLines(file.path(OUT,"LRG_final_symbols.txt"))

set.seed(20260922); N <- 3000
GS <- list("★ 乳酸化 LRG(332 底物蛋白)"=LRG, "乳酸化调控因子(writer/eraser)"=LACREG,
           "癫痫离子通道基因"=ION, "突触/囊泡基因"=SYN, "糖酵解通路"=GLYC)
res <- data.frame()
for(nm in names(GS)){
  hit <- toupper(tt$symbol) %in% toupper(unique(GS[[nm]]))
  k <- sum(hit); if(k<4) { cat(sprintf("  %-28s 检出=%d 跳过\n",nm,k)); next }
  obs <- sum(is_deg[hit]); nl <- replicate(N, sum(is_deg[sample.int(nrow(tt),k)]))
  res <- rbind(res, data.frame(gene_set=nm, n=k, n_sig=obs,
    Z=round((obs-mean(nl))/sd(nl),2), empP=signif((sum(nl>=obs)+1)/(N+1),3)))
}
cat("\n=== (a) bulk 转录组: 各基因集富集 ===\n")
print(res[order(-res$Z),], row.names=FALSE)

## ---------- (b) MR 命中率 ----------
cat("\n=== (b) MR: 各基因集命中率 vs 背景 ===\n")
lrg_ensg <- unique(sub("\\..*$","",trimws(readLines(file.path(OUT,"LRG_final_ensg.txt")))))
lrg_ensg <- lrg_ensg[nzchar(lrg_ensg)]
sym2ensg <- list()
for(s in c(ION,SYN,LACREG,GLYC)) sym2ensg[[s]] <- NA
# 用 ensg2symbol 反查
e2s <- sapply(names(map), function(g){ v<-map[[g]]
  if(is.null(v)) NA_character_ else if(is.list(v)) v$symbol[1] else NA_character_})
names(e2s) <- names(map)
mrres <- data.frame()
for(gk in c("Epilepsy_all","Focal_epilepsy","GGE")){
  f <- file.path(OUT, paste0("MRbg_gene_",gk,".csv"))
  if(!file.exists(f)) next
  gl <- read.csv(f, stringsAsFactors=FALSE)
  gl$symbol <- e2s[gl$ensg]
  bg_rate <- mean(gl$ratio_p < 0.05)
  cat(sprintf("\n[%s]  背景命中率=%.2f%%  (n=%d)\n", gk, 100*bg_rate, nrow(gl)))
  for(nm in names(GS)){
    genes <- unique(GS[[nm]])
    sub <- gl[toupper(gl$symbol) %in% toupper(genes),]
    if(nrow(sub)<4){ cat(sprintf("   %-28s 可检验=%d 跳过\n",nm,nrow(sub))); next }
    hit <- sum(sub$ratio_p<0.05); rate <- 100*hit/nrow(sub)
    others <- gl[!(toupper(gl$symbol) %in% toupper(genes)),]
    ft <- fisher.test(matrix(c(hit, nrow(sub)-hit, sum(others$ratio_p<0.05),
             nrow(others)-sum(others$ratio_p<0.05)), nrow=2, byrow=TRUE), alternative="greater")
    cat(sprintf("   %-28s %2d/%-4d = %5.1f%%   Fisher P=%.3f %s\n", nm, hit, nrow(sub), rate,
                ft$p.value, ifelse(ft$p.value<0.05,"-> 富集","" )))
    mrres <- rbind(mrres, data.frame(gwas=gk, gene_set=nm, n=nrow(sub), hit=hit,
                                     rate=round(rate,1), bg=round(100*bg_rate,1), fisher_p=signif(ft$p.value,3)))
  }
}
write.csv(res, file.path(OUT,"geneset_bulk_pathway_test.csv"), row.names=FALSE)
write.csv(mrres, file.path(OUT,"geneset_MR_pathway_test.csv"), row.names=FALSE)
cat("\n[保存] geneset_bulk_pathway_test.csv / geneset_MR_pathway_test.csv\n")
