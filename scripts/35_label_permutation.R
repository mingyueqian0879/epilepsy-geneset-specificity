# ============================================================
# Step 35: 校准检验 —— 标签置换 (label permutation)
#  H0: 病例/对照标签与表达无关
#  对每个数据集: 随机打乱标签 -> 跑完整 limma -> 记录 DEG 数
#              -> 同时对基因集做匹配置换检验 -> 记录 FDR<0.05 个数
#  若真实标签的结果落在零分布内 -> 该数据集的"发现"不可信
# ============================================================
suppressMessages({library(limma); library(edgeR)})
set.seed(20260922)
NLABEL <- 100        # 标签置换次数
NPERM  <- 500        # 基因集置换次数(每次标签置换内)

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
GS <- c(read_gmt(file.path(MS, "h.all.v2024.1.Hs.symbols.gmt")),
        read_gmt(file.path(MS, "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")))
cat("[gene sets]", length(GS), "\n")

mapdf <- read.csv(file.path(MOUT, "ensg2symbol_gencode.csv"), stringsAsFactors = FALSE)
ENS2SYM <- setNames(mapdf$symbol, mapdf$ensg)

matched_Z <- function(is_deg, ave, hit, nperm) {
  k <- sum(hit)
  if (k < 5 || sum(is_deg) < 3) return(NA_real_)
  obs <- sum(is_deg[hit])
  br <- unique(quantile(ave, probs = seq(0, 1, 0.05), na.rm = TRUE))
  if (length(br) < 3) return(NA_real_)
  st <- cut(ave, breaks = br, include.lowest = TRUE, labels = FALSE)
  st[is.na(st)] <- 1L
  need <- tabulate(st[hit], nbins = max(st))
  idx_by_st <- split(seq_along(st), st)
  sims <- numeric(nperm)
  for (i in seq_len(nperm)) {
    tot <- 0L
    for (s in seq_along(need)) {
      m <- need[s]; if (m == 0) next
      pool <- idx_by_st[[as.character(s)]]
      if (is.null(pool) || !length(pool)) next
      sel <- if (length(pool) <= m) pool else sample(pool, m)
      tot <- tot + sum(is_deg[sel])
    }
    sims[i] <- tot
  }
  if (sd(sims) == 0) return(NA_real_)
  (obs - mean(sims)) / sd(sims)
}

DATASETS <- c("GSE256068", "GSE213488", "GSE128300", "GSE272142")
out_rows <- list()

for (ds in DATASETS) {
  cat("\n", strrep("=", 76), "\n [", ds, "] 标签置换\n", sep = "")
  E <- as.matrix(read.delim(file.path(MOUT, paste0(ds, "_expr.tsv")),
                            row.names = 1, check.names = FALSE))
  dd <- read.csv(file.path(MOUT, paste0(ds, "_design.csv")), encoding = "UTF-8",
                 stringsAsFactors = FALSE)
  dd <- dd[dd$sample %in% colnames(E), ]
  E <- E[, dd$sample, drop = FALSE]
  is_counts <- max(E, na.rm = TRUE) > 100

  ## 预处理一次(标准化与标签无关的部分)
  if (is_counts) {
    keep <- rowSums(E >= 5) >= max(3, 0.2 * ncol(E))
    E <- E[keep, , drop = FALSE]
    dge <- DGEList(counts = round(E)); dge <- calcNormFactors(dge)
  } else {
    E <- E[rowSums(is.finite(E)) == ncol(E), , drop = FALSE]
  }

  run_once <- function(labels, do_gs = FALSE) {
    grp <- factor(labels, levels = c("Control", "Case"))
    design <- model.matrix(~ grp)
    if (is_counts) {
      vm <- voom(dge, design, plot = FALSE)
      fit <- eBayes(lmFit(vm, design), trend = TRUE, robust = TRUE)
    } else {
      fit <- eBayes(lmFit(E, design), trend = TRUE, robust = TRUE)
    }
    tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
    tt$symbol <- ENS2SYM[sub("\\..*$", "", rownames(tt))]
    tt <- tt[!is.na(tt$symbol), ]
    is_deg <- tt$adj.P.Val < 0.05
    if (!do_gs) return(sum(is_deg))
    zs <- sapply(GS, function(g) matched_Z(is_deg, tt$AveExpr,
                                           toupper(tt$symbol) %in% toupper(g), NPERM))
    zs <- zs[!is.na(zs)]
    if (!length(zs)) return(c(n_sig = 0, n_test = 0))
    p <- 2 * pnorm(-abs(zs))
    fdr <- p.adjust(p, "BH")
    c(n_sig = sum(fdr < 0.05), n_test = length(zs))
  }

  ## 真实标签
  obs_deg <- run_once(dd$group, do_gs = FALSE)
  obs_gs <- run_once(dd$group, do_gs = TRUE)
  cat(sprintf("  真实标签: DEG=%d (%.1f%%)  基因集FDR<0.05=%d/%d\n",
              obs_deg, 100 * obs_deg / nrow(E), obs_gs["n_sig"], obs_gs["n_test"]))

  ## 标签置换零分布
  null_deg <- numeric(NLABEL); null_gs <- numeric(NLABEL)
  for (i in seq_len(NLABEL)) {
    lab <- sample(dd$group)
    null_deg[i] <- run_once(lab, do_gs = FALSE)
    r <- run_once(lab, do_gs = TRUE)
    null_gs[i] <- r["n_sig"]
    if (i %% 20 == 0) cat(sprintf("     label perm %d/%d ...\n", i, NLABEL))
  }

  p_deg <- (sum(null_deg >= obs_deg) + 1) / (NLABEL + 1)
  p_gs  <- (sum(null_gs >= obs_gs["n_sig"]) + 1) / (NLABEL + 1)
  cat(sprintf("  零分布 DEG: 中位=%.0f  95%%=%.0f  最大=%.0f\n",
              median(null_deg), quantile(null_deg, .95), max(null_deg)))
  cat(sprintf("  >>> DEG 经验P = %.4f        (真实DEG=%d vs 零中位=%.0f)\n",
              p_deg, obs_deg, median(null_deg)))
  cat(sprintf("  >>> 基因集 经验P = %.4f     (真实FDR<0.05=%d vs 零中位=%.0f)\n",
              p_gs, obs_gs["n_sig"], median(null_gs)))

  out_rows[[ds]] <- data.frame(
    dataset = ds, n_case = sum(dd$group == "Case"), n_ctrl = sum(dd$group == "Control"),
    n_genes = nrow(E), obs_deg = obs_deg,
    obs_deg_pct = round(100 * obs_deg / nrow(E), 2),
    null_deg_median = median(null_deg), null_deg_p95 = quantile(null_deg, .95),
    p_deg = p_deg, obs_gs_sig = obs_gs["n_sig"], gs_tested = obs_gs["n_test"],
    null_gs_median = median(null_gs), p_gs = p_gs)
  saveRDS(list(null_deg = null_deg, null_gs = null_gs, obs_deg = obs_deg,
               obs_gs = obs_gs), file.path(OUT, paste0("labelperm_", ds, ".rds")))
}

R <- do.call(rbind, out_rows)
write.csv(R, file.path(OUT, "label_permutation.csv"), row.names = FALSE)
cat("\n", strrep("=", 78), "\n 标签置换校准汇总\n", strrep("=", 78), "\n", sep = "")
print(R, row.names = FALSE)
cat("\n[saved] label_permutation.csv\n")
