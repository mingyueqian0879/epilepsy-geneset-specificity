# ============================================================
# Step 35b: 标签置换校准 (快速版) —— 只检验 DEG 数
#  H0: 分组标签与表达无关
#  对每个数据集打乱标签 N 次 -> 跑完整 limma -> 记录 DEG 数
#  与真实标签比较 -> 经验 P
# ============================================================
suppressMessages({library(limma); library(edgeR)})
set.seed(20260922)
NLABEL <- 200

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
MOUT <- file.path(BASE, "out", "multi")
OUT  <- file.path(BASE, "out")

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
  cat("\n", strrep("=", 74), "\n [", ds, "]  genes=", nrow(E), " samples=", ncol(E), "\n", sep = "")

  n_deg_for <- function(labels) {
    grp <- factor(labels, levels = c("Control", "Case"))
    design <- model.matrix(~ grp)
    if (is_counts) {
      vm <- voom(dge, design, plot = FALSE)
      fit <- eBayes(lmFit(vm, design), trend = TRUE, robust = TRUE)
    } else {
      fit <- eBayes(lmFit(E, design), trend = TRUE, robust = TRUE)
    }
    tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
    sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  }

  obs <- n_deg_for(dd$group)
  null <- numeric(NLABEL)
  for (i in seq_len(NLABEL)) {
    null[i] <- n_deg_for(sample(dd$group))
    if (i %% 50 == 0) cat(sprintf("   perm %d/%d  当前零中位=%.0f\n", i, NLABEL, median(null[1:i])))
  }
  p <- (sum(null >= obs) + 1) / (NLABEL + 1)
  cat(sprintf("  真实 DEG = %d (%.1f%%)\n", obs, 100 * obs / nrow(E)))
  cat(sprintf("  零分布   : 中位=%.0f  IQR=%.0f-%.0f  95%%=%.0f  最大=%.0f\n",
              median(null), quantile(null, .25), quantile(null, .75),
              quantile(null, .95), max(null)))
  cat(sprintf("  零分布 DEG 率中位 = %.2f%%\n", 100 * median(null) / nrow(E)))
  cat(sprintf("  >>> 经验 P = %.4f   %s\n", p,
              ifelse(p < 0.05, "显著(DE信号真实存在)", "不显著(未检出信号)")))
  cat(sprintf("  >>> 倍数 = %.1fx 零中位\n", obs / max(median(null), 1)))

  rows[[ds]] <- data.frame(dataset = ds, n_genes = nrow(E),
    n_case = sum(dd$group == "Case"), n_ctrl = sum(dd$group == "Control"),
    obs_deg = obs, obs_pct = round(100 * obs / nrow(E), 2),
    null_median = median(null), null_median_pct = round(100 * median(null) / nrow(E), 2),
    null_p95 = quantile(null, .95), fold = round(obs / max(median(null), 1), 2),
    emp_P = p)
  saveRDS(list(null = null, obs = obs), file.path(OUT, paste0("labelperm_", ds, ".rds")))
}

R <- do.call(rbind, rows)
write.csv(R, file.path(OUT, "label_permutation_200perm.csv"), row.names = FALSE)
cat("\n", strrep("=", 84), "\n 标签置换校准汇总\n", strrep("=", 84), "\n", sep = "")
print(R, row.names = FALSE)
cat("\n[saved] label_permutation_200perm.csv\n")

