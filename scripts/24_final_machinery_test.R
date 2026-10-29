# ============================================================
# Step ②-F: 最终验证 —— "乳酸化机器"假说的稳健性
#  核心发现: 332 底物蛋白列表无信号, 但 writer/eraser 有信号
#  需验证: (1) 扩展基因集后是否保持
#          (2) 多重检验校正后是否成立
#          (3) 是否被细胞组成混杂解释
#          (4) LOO 敏感性
# ============================================================
suppressMessages({library(limma); library(jsonlite)})
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); OUT <- file.path(BASE,"out")
set.seed(20260922); N <- 5000

expr_df <- read.delim(file.path(OUT,"expr_matrix.tsv"), row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
expr <- as.matrix(expr_df); storage.mode(expr) <- "double"
meta <- read.csv(file.path(OUT,"design_aligned.csv"), stringsAsFactors=FALSE, fileEncoding="UTF-8")
rownames(meta) <- meta$matrix_col; meta <- meta[colnames(expr),]
map <- fromJSON(file.path(OUT,"ensg2symbol.json"))
sym <- sapply(sub("\\..*$","",rownames(expr)), function(g){ v<-map[[g]]
  if(is.null(v)) NA_character_ else if(is.list(v)) v$symbol[1] else NA_character_})

## ---------- 基因集定义 ----------
WRITER <- c("EP300","CREBBP","KAT2A","KAT2B","KAT5","KAT6A","KAT6B","KAT7","KAT8","ELP3","NAA10","KAT14","KAT2B")
ERASER <- c("HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10","HDAC11",
            "SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7")
READER <- c("BRD2","BRD3","BRD4","BRDT","YEATS2","ZZZ3","MLLT3","YEATS4","SP140")
MACHINERY <- unique(c(WRITER, ERASER, READER))
LACTMETAB <- c("LDHA","LDHB","LDHC","SLC16A1","SLC16A3","SLC16A7","HIF1A","PDK1","PKM","MPC1","MPC2","BSG")
SUBSTRATE <- readLines(file.path(OUT,"LRG_final_symbols.txt"))

cat(sprintf("[基因集规模] writer=%d eraser=%d reader=%d -> machinery=%d ; 底物=%d\n",
            length(unique(WRITER)),length(unique(ERASER)),length(unique(READER)),
            length(MACHINERY),length(SUBSTRATE)))

## ---------- 模型 ----------
keep <- meta$disease %in% c("TLE-HS","ControlCortex")
e <- expr[,keep]; m <- meta[keep,]
m$grp <- factor(ifelse(m$disease=="TLE-HS","Case","Ctrl"), levels=c("Ctrl","Case"))
z <- t(scale(t(e)))
MARK <- list(astro=c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","GJA1","S100B","ALDOC","GLUL","SOX9","NFIA","NFIB","CD44","VIM"),
 micro=c("TMEM119","P2RY12","CX3CR1","AIF1","CD68","CSF1R","TREM2","TYROBP","C1QA","C1QB","C1QC","ITGAM","SPI1","CD14"),
 oligo=c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","SOX10","OLIG1","OLIG2"),
 neuro=c("RBFOX3","SYP","SNAP25","NEFL","NEFM","MAP2","TUBB3","GRIN1","SLC17A7","GAD1","GAD2","SYT1","NEFH","STMN2","UCHL1"),
 immune=c("CD3D","CD3E","CD2","CD8A","CD4","MS4A1","CD79A","CD79B","IGHM","CD37","SELPLG","IKZF1","WAS","MNDA","PYCARD"))
for(nm in names(MARK)){ idx<-which(toupper(sym) %in% MARK[[nm]]); m[[paste0("sc_",nm)]]<-colMeans(z[idx,,drop=FALSE],na.rm=TRUE) }

fitmod <- function(form){
  design <- model.matrix(form, data=m)
  fit <- eBayes(lmFit(e,design), trend=TRUE, robust=TRUE)
  tt <- topTable(fit, coef="grpCase", number=Inf, sort.by="P")
  tt$symbol <- sym[rownames(tt)]          # 关键: 必须带基因名
  tt
}
ttA <- fitmod(~ grp + factor(tissue))
ttB <- fitmod(~ grp + factor(tissue) + sc_astro + sc_micro + sc_oligo + sc_neuro + sc_immune)

enrich <- function(tt, genes, lab){
  is_deg <- tt$adj.P.Val < 0.05
  hit <- toupper(tt$symbol) %in% toupper(unique(genes))
  k <- sum(hit); if(k<4) return(NULL)
  obs <- sum(is_deg[hit])
  nl <- replicate(N, sum(is_deg[sample.int(nrow(tt),k)]))
  data.frame(set=lab, n=k, n_sig=obs, pct=round(100*obs/k,1),
             Z=round((obs-mean(nl))/sd(nl),2), empP=signif((sum(nl>=obs)+1)/(N+1),3))
}
SETS <- list("★ 乳酸化机器(writer+eraser+reader)"=MACHINERY, "writer(乙酰转移酶)"=WRITER,
             "eraser(HDAC+SIRT)"=ERASER, "reader(BRD/YEATS)"=READER,
             "乳酸代谢"=LACTMETAB, "★ 底物蛋白列表(332)"=SUBSTRATE)

cat("\n========== 模型A: ~grp + tissue (未校正) ==========\n")
la <- lapply(names(SETS), function(n) enrich(ttA, SETS[[n]], n))
la <- la[!vapply(la, is.null, logical(1))]
ra <- do.call(rbind, la)
ra$bh <- p.adjust(ra$empP, "BH")
print(ra[order(-ra$Z),], row.names=FALSE)
cat(sprintf("  [全基因组 DEG 基线 = %.1f%%]\n", 100*mean(ttA$adj.P.Val<0.05)))

cat("\n========== 模型B: 细胞组成校正后 ==========\n")
lb <- lapply(names(SETS), function(n) enrich(ttB, SETS[[n]], n))
lb <- lb[!vapply(lb, is.null, logical(1))]
rb <- do.call(rbind, lb)
rb$bh <- p.adjust(rb$empP, "BH")
print(rb[order(-rb$Z),], row.names=FALSE)
cat(sprintf("  [全基因组 DEG 基线 = %.1f%%]\n", 100*mean(ttB$adj.P.Val<0.05)))

## ---------- LOO 敏感性 (乳酸化机器) ----------
cat("\n========== LOO 敏感性: 逐个剔除基因(模型A) ==========\n")
gene_loo <- sapply(MACHINERY, function(g){
  s <- setdiff(MACHINERY, g)
  r <- enrich(ttA, s, g)
  if(is.null(r)) return(c(NA,NA)) else c(r$Z, r$empP)
})
loo <- data.frame(dropped=colnames(gene_loo), Z=round(gene_loo[1,],2), P=signif(gene_loo[2,],3))
loo <- loo[order(-loo$Z),]
cat(sprintf("  LOO Z 范围: %.2f ~ %.2f   中位 %.2f\n", min(loo$Z,na.rm=TRUE), max(loo$Z,na.rm=TRUE), median(loo$Z,na.rm=TRUE)))
print(head(loo,8), row.names=FALSE)
print(tail(loo,5), row.names=FALSE)

## ---------- 逐基因效应 (乳酸化机器) ----------
cat("\n========== 乳酸化机器基因的逐个效应 (模型A) ==========\n")
mm <- ttA[toupper(ttA$symbol) %in% toupper(MACHINERY), c("symbol","logFC","AveExpr","P.Value","adj.P.Val")]
mm$cls <- ifelse(mm$symbol %in% WRITER,"writer", ifelse(mm$symbol %in% ERASER,"eraser","reader"))
print(mm[order(mm$P.Value),], row.names=FALSE, digits=3)

write.csv(ra, file.path(OUT,"final_lactylation_machinery_modelA.csv"), row.names=FALSE)
write.csv(rb, file.path(OUT,"final_lactylation_machinery_modelB.csv"), row.names=FALSE)
write.csv(mm, file.path(OUT,"final_machinery_genes.csv"), row.names=FALSE)
cat("\n[保存] final_lactylation_machinery_*.csv  final_machinery_genes.csv\n")
