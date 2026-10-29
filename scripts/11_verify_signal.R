# ============================================================
# Step 5: 关键验证 —— 是真信号还是混杂?
#  1) 把显著 LRG 映射成基因名，看是什么
#  2) 阴性对照：1000 个随机基因集(同等大小) -> 零分布
#  3) 替代对比：TLE-HS vs 其他癫痫病理（都是手术标本，控制取材混杂）
# ============================================================
suppressMessages({library(limma); library(jsonlite)})
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); OUT <- file.path(BASE, "out")

expr_df <- read.delim(file.path(OUT, "expr_matrix.tsv"), row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
expr <- as.matrix(expr_df); storage.mode(expr) <- "double"
meta <- read.csv(file.path(OUT, "design_aligned.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
rownames(meta) <- meta$matrix_col; meta <- meta[colnames(expr), ]

## 基因名映射
map <- fromJSON(file.path(OUT, "ensg2symbol.json"))
sym <- sapply(sub("\\..*$", "", rownames(expr)), function(g) {
  v <- map[[g]]; if (is.null(v)) NA_character_ else if (is.list(v)) v$symbol[1] else NA_character_
})
names(sym) <- rownames(expr)

lrg_ensg <- unique(sub("\\..*$", "", trimws(readLines(file.path(OUT, "LRG_final_ensg.txt")))))
lrg_ensg <- lrg_ensg[nzchar(lrg_ensg)]
cat("[LRG] ", length(lrg_ensg), " 个 ENSG\n")

## ---------- 1. 主对比 + 显著 LRG 身份 ----------
run_limma <- function(case_lab, ctrl_lab, tag) {
  keep <- meta$disease %in% c(case_lab, ctrl_lab)
  e <- expr[, keep]; m <- meta[keep, ]
  m$grp <- factor(ifelse(m$disease == case_lab, "Case", "Ctrl"), levels = c("Ctrl", "Case"))
  m$tis <- factor(m$tissue)
  design <- model.matrix(~ grp + tis, data = m)
  fit <- eBayes(lmFit(e, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
  tt$ensg <- sub("\\..*$", "", rownames(tt)); tt$symbol <- sym[rownames(tt)]
  tt$is_lrg <- tt$ensg %in% lrg_ensg
  cat(sprintf("\n[%s]  Case=%d Ctrl=%d\n", tag, sum(m$grp=="Case"), sum(m$grp=="Ctrl")))
  tt
}

tt <- run_limma("TLE-HS", "ControlCortex", "TLE-HS vs ControlCortex")
deg <- sum(tt$adj.P.Val < 0.05)
cat(sprintf("  总基因=%d  DEG(adj.P<0.05)=%d (%.1f%%)\n", nrow(tt), deg, 100*deg/nrow(tt)))

cat("\n--- 全部基因中 top15 (看是什么在驱动差异) ---\n")
print(head(tt[, c("symbol","logFC","AveExpr","adj.P.Val")], 15), row.names = FALSE, digits = 3)

sig_lrg <- tt[tt$is_lrg & tt$adj.P.Val < 0.05, ]
cat(sprintf("\n--- 显著 LRG: %d 个, top20 ---\n", nrow(sig_lrg)))
print(head(sig_lrg[, c("symbol","logFC","AveExpr","adj.P.Val")], 20), row.names = FALSE, digits = 3)

## ---------- 2. 阴性对照：随机基因集零分布 ----------
cat("\n======================================================\n")
cat(" 阴性对照: 1000 个随机基因集 (n=LRG检出数)\n")
cat("======================================================\n")
lrg_present <- tt$ensg[tt$is_lrg]
n_lrg <- length(lrg_present)
cat(sprintf("  LRG 可检验数 = %d\n", n_lrg))

set.seed(20260922)
N <- 1000
rand_counts <- numeric(N)
idx_all <- seq_len(nrow(tt))
is_deg <- tt$adj.P.Val < 0.05
for (i in seq_len(N)) {
  s <- sample(idx_all, n_lrg)
  rand_counts[i] <- sum(is_deg[s])
}
obs <- sum(is_deg[tt$is_lrg])
cat(sprintf("\n  观测显著 LRG 数        = %d\n", obs))
cat(sprintf("  随机基因集显著数 均值  = %.1f  (SD %.1f)\n", mean(rand_counts), sd(rand_counts)))
cat(sprintf("  随机分布 2.5%%~97.5%%   = %.0f ~ %.0f\n",
            quantile(rand_counts, .025), quantile(rand_counts, .975)))
cat(sprintf("  经验 P (观测 > 随机)   = %.4f\n", (sum(rand_counts >= obs) + 1) / (N + 1)))
cat(sprintf("  Z 分数                 = %.3f\n", (obs - mean(rand_counts)) / sd(rand_counts)))

## ---------- 3. 替代对比: 手术标本之间 (控制取材混杂) ----------
cat("\n======================================================\n")
cat(" 替代对比: 手术标本之间 (控制'手术 vs 尸检'混杂)\n")
cat("======================================================\n")
meta$grp2 <- ifelse(meta$disease == "TLE-HS", "TLE-HS", "OtherEpi")
for (cc in c("FCD2b", "TSC", "FCD2a")) {
  tt2 <- run_limma("TLE-HS", cc, paste0("TLE-HS vs ", cc))
  d2 <- sum(tt2$adj.P.Val < 0.05)
  o2 <- sum(tt2$is_lrg & tt2$adj.P.Val < 0.05)
  n2 <- sum(tt2$is_lrg)
  a <- o2; b <- n2 - o2
  c_ <- sum(!tt2$is_lrg & tt2$adj.P.Val < 0.05); d_ <- sum(!tt2$is_lrg & tt2$adj.P.Val >= 0.05)
  ft <- fisher.test(matrix(c(a,b,c_,d_), nrow=2, byrow=TRUE), alternative="greater")
  cat(sprintf("  DEG=%d (%.1f%%)  LRG检验=%d  显著LRG=%d  Fisher OR=%.2f P=%.3g\n",
              d2, 100*d2/nrow(tt2), n2, o2, ft$estimate, ft$p.value))
}

write.csv(tt, file.path(OUT, "limma_main_annotated.csv"), row.names = FALSE)
write.csv(sig_lrg, file.path(OUT, "significant_LRGs.csv"), row.names = FALSE)
saveRDS(list(rand=rand_counts, obs=obs), file.path(OUT, "negative_control.rds"))
cat("\n[保存] limma_main_annotated.csv / significant_LRGs.csv / negative_control.rds\n")
