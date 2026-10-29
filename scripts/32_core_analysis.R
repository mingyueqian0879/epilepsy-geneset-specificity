# ============================================================
# Step 32: 核心分析 —— 多数据集 x 多基因集 特异性检验矩阵
#
#  每个数据集:
#   1) limma 差异表达 (counts -> voom; log2TPM -> 直接)
#   2) 全基因组 DEG 率
#   3) 系统性基因集检验: Hallmark(50) + KEGG + 热点基因集 + 阳性对照
#   4) 零分布: 表达量分层匹配置换 (控制"高表达基因更易显著"偏倚)
#   5) 输出 Z 分数矩阵
# ============================================================
suppressMessages({library(limma); library(edgeR)})
set.seed(20260922)
NPERM <- 2000

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
MOUT <- file.path(BASE, "out", "multi")
OUT  <- file.path(BASE, "out")

read_gmt <- function(f) {
  out <- list()
  for (l in readLines(f)) {
    p <- strsplit(l, "\t")[[1]]
    if (length(p) >= 3) out[[p[1]]] <- unique(p[-c(1, 2)])
  }
  out
}
MS <- file.path(BASE, "data", "msigdb")
GS <- list(HALLMARK = read_gmt(file.path(MS, "h.all.v2024.1.Hs.symbols.gmt")),
           KEGG     = read_gmt(file.path(MS, "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")))

HOT <- list(
 "Lactylation_LRG332" = readLines(file.path(OUT, "LRG_final_symbols.txt")),
 "Ferroptosis_core" = c("GPX4","ACSL4","SLC7A11","FTH1","FTL","NCOA4","TFRC","ALOX15","SAT1",
   "GCLC","GCLM","SLC3A2","PCBP1","PCBP2","HMOX1","NFE2L2","KEAP1","ATF4","DDIT3","CHAC1",
   "SLC40A1","STEAP3","CISD1","FANCD2","PROM2","ABCC1","LPCAT3","ALOX12","ALOX5","POR",
   "CYBB","NOS2","NOX1","NOX4","VDAC2","VDAC3","G6PD","SLC1A5","GLS2","PHKG2"),
 "Cuproptosis_core" = c("FDX1","LIAS","LIPT1","DLD","DLAT","PDHA1","PDHB","MTF1","GLS",
   "CDKN2A","SLC31A1","ATP7A","ATP7B","DLST","GCSH","PDHX","LIPT2","TIMM8B","NDUFA13","SDHB"),
 "Disulfidptosis_core" = c("SLC7A11","SLC3A2","SLC7A9","GYS1","RPN1","NCKAP1","OXSM","NCKAP1L",
   "ACTN4","ACTB","ACTN1","ACTN2","MYH9","MYH10","MYL6","FLNA","FLNB","TLN1","CAPZB","CD2AP",
   "INF2","DSTN","PGM1","G6PD","TKT","RPN2","SLC7A5","SLC7A8"),
 "Pyroptosis_core" = c("NLRP3","CASP1","CASP4","CASP5","CASP8","GSDMD","GSDME","IL1B","IL18",
   "AIM2","NLRC4","PYCARD","NLRP1","NLRP6","NLRP9","NLRC5","MEFV","GSDMA","GSDMB","GSDMC"),
 "m6A_regulators" = c("METTL3","METTL14","METTL16","WTAP","VIRMA","ZC3H13","RBM15","RBM15B",
   "YTHDF1","YTHDF2","YTHDF3","YTHDC1","YTHDC2","FTO","ALKBH5","HNRNPA2B1","HNRNPC",
   "IGF2BP1","IGF2BP2","IGF2BP3","ELAVL1","LRPPRC","FMR1","PRRC2A"),
 "Autophagy_core" = c("ATG3","ATG4A","ATG4B","ATG4C","ATG4D","ATG5","ATG7","ATG9A","ATG9B",
   "ATG10","ATG12","ATG13","ATG14","ATG16L1","ATG16L2","BECN1","MAP1LC3A","MAP1LC3B",
   "MAP1LC3C","SQSTM1","ULK1","ULK2","ULK3","ULK4","PIK3C3","LAMP1","LAMP2","TFEB","MTOR","AMBRA1"),
 "Glycolysis_core" = c("HK1","HK2","HK3","GPI","PFKL","PFKM","PFKP","ALDOA","ALDOB","ALDOC",
   "TPI1","GAPDH","PGK1","PGAM1","PGAM2","ENO1","ENO2","ENO3","PKM","PKLR","LDHA","LDHB","LDHC"),
 "POS_Astrocyte"  = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","GJA1","S100B","ALDOC","GLUL",
                      "SOX9","NFIA","NFIB","CD44","VIM"),
 "POS_Microglia"  = c("TMEM119","P2RY12","CX3CR1","AIF1","CD68","CSF1R","TREM2","TYROBP",
                      "C1QA","C1QB","C1QC","ITGAM","SPI1","CD14"),
 "POS_Neuron"     = c("RBFOX3","SYP","SNAP25","NEFL","NEFM","MAP2","TUBB3","GRIN1","SLC17A7",
                      "GAD1","GAD2","SYT1","NEFH","STMN2","UCHL1"),
 "POS_Oligodendrocyte" = c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","SOX10","OLIG1","OLIG2"))
GS[["HOT"]] <- HOT
cat(sprintf("[gene sets] HALLMARK=%d KEGG=%d HOT=%d total=%d\n",
            length(GS$HALLMARK), length(GS$KEGG), length(GS$HOT),
            sum(sapply(GS, length))))

mapdf <- read.csv(file.path(MOUT, "ensg2symbol_gencode.csv"), stringsAsFactors = FALSE)
ENS2SYM <- setNames(mapdf$symbol, mapdf$ensg)

## 表达量分层匹配置换: 对每个基因集返回 Z 与经验 P
matched_test <- function(is_deg, ave, hit, nperm = NPERM) {
  k <- sum(hit)
  if (k < 5) return(c(n = k, obs = NA, exp = NA, Z = NA, P = NA))
  obs <- sum(is_deg[hit])
  ## 按表达量分 20 层
  br <- unique(quantile(ave, probs = seq(0, 1, 0.05), na.rm = TRUE))
  if (length(br) < 3) br <- unique(quantile(ave, probs = seq(0, 1, 0.5), na.rm = TRUE))
  st <- cut(ave, breaks = br, include.lowest = TRUE, labels = FALSE)
  st[is.na(st)] <- 1L
  ## hit 基因在各层的数量
  need <- tabulate(st[hit], nbins = max(st))
  idx_by_st <- split(seq_along(st), st)
  sims <- numeric(nperm)
  for (i in seq_len(nperm)) {
    tot <- 0L
    for (s in seq_along(need)) {
      m <- need[s]
      if (m == 0) next
      pool <- idx_by_st[[as.character(s)]]
      if (is.null(pool) || length(pool) == 0) next
      sel <- if (length(pool) <= m) pool else sample(pool, m)
      tot <- tot + sum(is_deg[sel])
    }
    sims[i] <- tot
  }
  c(n = k, obs = obs, exp = mean(sims), Z = (obs - mean(sims)) / sd(sims),
    P = (sum(sims >= obs) + 1) / (nperm + 1))
}

DATASETS <- c("GSE256068", "GSE213488", "GSE128300", "GSE272142")
all_rows <- list()
summary_ds <- data.frame()

for (ds in DATASETS) {
  cat("\n", strrep("=", 78), "\n [", ds, "]\n", sep = "")
  e <- as.matrix(read.delim(file.path(MOUT, paste0(ds, "_expr.tsv")),
                            row.names = 1, check.names = FALSE))
  dd <- read.csv(file.path(MOUT, paste0(ds, "_design.csv")), encoding = "UTF-8",
                 stringsAsFactors = FALSE)
  dd <- dd[dd$sample %in% colnames(e), ]
  e <- e[, dd$sample, drop = FALSE]

  ## 去低表达
  if (max(e, na.rm = TRUE) > 100) {          # 判定为 counts
    keep <- rowSums(e >= 5) >= max(3, 0.2 * ncol(e))
    e <- e[keep, , drop = FALSE]
    dge <- DGEList(counts = round(e))
    dge <- calcNormFactors(dge)
    grp <- factor(dd$group, levels = c("Control", "Case"))
    design <- model.matrix(~ grp)
    vm <- voom(dge, design, plot = FALSE)
    fit <- eBayes(lmFit(vm, design), trend = TRUE, robust = TRUE)
    cat("   [mode] counts -> voom\n")
  } else {
    e <- e[rowSums(is.finite(e)) == ncol(e), , drop = FALSE]
    grp <- factor(dd$group, levels = c("Control", "Case"))
    design <- model.matrix(~ grp)
    fit <- eBayes(lmFit(e, design), trend = TRUE, robust = TRUE)
    cat("   [mode] log2 TPM -> limma\n")
  }
  tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
  tt$ensg <- sub("\\..*$", "", rownames(tt))
  tt$symbol <- ENS2SYM[tt$ensg]
  tt <- tt[!is.na(tt$symbol), ]

  is_deg <- tt$adj.P.Val < 0.05
  deg_rate <- mean(is_deg)
  cat(sprintf("   genes tested = %d   DEG(adj.P<0.05) = %d  (%.1f%%)\n",
              nrow(tt), sum(is_deg), 100 * deg_rate))
  summary_ds <- rbind(summary_ds, data.frame(
    dataset = ds, n_case = sum(dd$group == "Case"), n_ctrl = sum(dd$group == "Control"),
    n_genes = nrow(tt), n_deg = sum(is_deg), deg_rate = round(100 * deg_rate, 2)))

  for (fam in names(GS)) {
    for (gn in names(GS[[fam]])) {
      g <- toupper(unique(GS[[fam]][[gn]]))
      hit <- toupper(tt$symbol) %in% g
      r <- matched_test(is_deg, tt$AveExpr, hit)
      all_rows[[length(all_rows) + 1]] <- data.frame(
        dataset = ds, family = fam, gene_set = gn,
        n = r["n"], obs = r["obs"], exp = round(r["exp"], 2),
        Z = round(r["Z"], 3), P = signif(r["P"], 3))
    }
  }
  cat(sprintf("   gene sets tested = %d\n", length(all_rows)))
}

R <- do.call(rbind, all_rows)
R$FDR <- NA
for (ds in DATASETS) {
  i <- R$dataset == ds & !is.na(R$P)
  R$FDR[i] <- p.adjust(R$P[i], "BH")
}
write.csv(R, file.path(OUT, "specificity_matrix.csv"), row.names = FALSE)
write.csv(summary_ds, file.path(OUT, "dataset_summary.csv"), row.names = FALSE)

cat("\n", strrep("=", 78), "\n 数据集汇总\n", strrep("=", 78), "\n", sep = "")
print(summary_ds, row.names = FALSE)

cat("\n", strrep("=", 78), "\n 各数据集: 基因集 Z 分数分布\n", strrep("=", 78), "\n", sep = "")
for (ds in DATASETS) {
  s <- R[R$dataset == ds & !is.na(R$Z), ]
  sig <- sum(s$FDR < 0.05, na.rm = TRUE)
  cat(sprintf("  %-12s  检验=%-4d  Z中位=%+.3f  Z均值=%+.3f  |Z|>2=%-3d  FDR<0.05=%-3d (%.1f%%)\n",
              ds, nrow(s), median(s$Z), mean(s$Z), sum(abs(s$Z) > 2), sig,
              100 * sig / nrow(s)))
}

cat("\n--- 热点基因集结果 (全部数据集) ---\n")
h <- R[R$family == "HOT", c("dataset", "gene_set", "n", "obs", "exp", "Z", "P", "FDR")]
print(h[order(h$gene_set, h$dataset), ], row.names = FALSE)
cat("\n[saved] specificity_matrix.csv / dataset_summary.csv\n")
