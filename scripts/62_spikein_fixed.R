# ============================================================
# Step 62: Spiked-in 模拟 (修正版)
#  修正: 经验 P 受置换次数限制(下限 1/(n+1)), 无法达到 BH 门槛。
#        改用置换零分布估计 Z, 再由正态近似计算 P (无分辨率限制)。
#  目的: 给出真实的检出率曲线, 验证 MDE 预测
# ============================================================
suppressMessages({library(limma)})
set.seed(20260922)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
MOUT <- file.path(BASE, "out", "multi")
OUT  <- file.path(BASE, "out")

E <- as.matrix(read.delim(file.path(MOUT, "GSE256068_expr.tsv"),
                          row.names = 1, check.names = FALSE))
dd <- read.csv(file.path(MOUT, "GSE256068_design.csv"), encoding = "UTF-8",
               stringsAsFactors = FALSE)
dd <- dd[dd$sample %in% colnames(E), ]
E <- E[, dd$sample, drop = FALSE]
grp <- factor(dd$group, levels = c("Control", "Case"))
cat(sprintf("[data] %d genes x %d samples (Case=%d, Control=%d)\n",
            nrow(E), ncol(E), sum(grp == "Case"), sum(grp == "Control")))

SETSIZE <- 100; NREP <- 30; NPERM <- 2000; M_TESTS <- 248
z_thresh <- qnorm(1 - (0.05 / M_TESTS) / 2)      # BH 秩1 对应的 Z 门槛

cat(sprintf("[thresholds] BH 秩1: P < %.2e  ->  Z >= %.2f\n",
            0.05 / M_TESTS, z_thresh))
cat(sprintf("[note] 2000 次置换的经验 P 下限 = %.5f, 无法达到上述门槛;\n", 1 / (NPERM + 1)))
cat("       故以置换零分布估计 Z, 再用正态近似求 P。\n\n")

matched_Z <- function(is_deg, ave, hit, nperm = NPERM) {
  k <- sum(hit); if (k < 5) return(NA_real_)
  obs <- sum(is_deg[hit])
  br <- unique(quantile(ave, probs = seq(0, 1, 0.05), na.rm = TRUE))
  st <- cut(ave, breaks = br, include.lowest = TRUE, labels = FALSE)
  st[is.na(st)] <- 1L
  need <- tabulate(st[hit], nbins = max(st))
  ib <- split(seq_along(st), st)
  sims <- replicate(nperm, {
    t <- 0L
    for (s in seq_along(need)) {
      m <- need[s]; if (!m) next
      pool <- ib[[as.character(s)]]
      if (is.null(pool) || !length(pool)) next
      t <- t + sum(is_deg[if (length(pool) <= m) pool else sample(pool, m)])
    }; t })
  if (sd(sims) == 0) return(NA_real_)
  (obs - mean(sims)) / sd(sims)
}

run_de <- function(Emat) {
  design <- model.matrix(~ grp)
  fit <- eBayes(lmFit(Emat, design), trend = TRUE, robust = TRUE)
  topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
}

DELTAS <- c(0, 0.15, 0.25, 0.35, 0.5, 0.75, 1.0)
res <- list()
for (dl in DELTAS) {
  cat(sprintf("[spike delta=%.2f] %d reps ...\n", dl, NREP))
  det <- 0; zs <- c(); exc <- c(); ob <- c(); ex <- c()
  for (r in seq_len(NREP)) {
    set.seed(20260922 + r)
    spiked <- sample(rownames(E), SETSIZE)
    Es <- E
    if (dl > 0) {
      ci <- which(grp == "Case")
      Es[spiked, ci] <- Es[spiked, ci] + dl
    }
    tt <- run_de(Es)
    isd <- tt$adj.P.Val < 0.05
    hit <- rownames(tt) %in% spiked
    k <- sum(hit); obs <- sum(isd[hit])
    # 零分布
    br <- unique(quantile(tt$AveExpr, probs = seq(0, 1, 0.05), na.rm = TRUE))
    st <- cut(tt$AveExpr, breaks = br, include.lowest = TRUE, labels = FALSE)
    st[is.na(st)] <- 1L
    need <- tabulate(st[hit], nbins = max(st)); ib <- split(seq_along(st), st)
    sims <- replicate(NPERM, {
      t <- 0L
      for (s in seq_along(need)) {
        m <- need[s]; if (!m) next
        pool <- ib[[as.character(s)]]
        if (is.null(pool) || !length(pool)) next
        t <- t + sum(isd[if (length(pool) <= m) pool else sample(pool, m)])
      }; t })
    if (sd(sims) > 0) {
      z <- (obs - mean(sims)) / sd(sims)
      zs <- c(zs, z); ob <- c(ob, obs); ex <- c(ex, mean(sims))
      exc <- c(exc, (obs - mean(sims)) / max(mean(sims), 0.5))
      if (z >= z_thresh) det <- det + 1
    }
  }
  res[[as.character(dl)]] <- data.frame(
    delta = dl, n = length(zs),
    detection_rate = round(100 * det / max(length(zs), 1), 1),
    median_Z = round(median(zs), 2),
    median_obs = round(median(ob), 1), median_exp = round(median(ex), 1),
    median_excess_pct = round(100 * median(exc), 1),
    fp_at_alpha05 = round(100 * mean(2 * pnorm(-abs(zs)) < 0.05), 1))
  cat(sprintf("   -> 检出率 %.1f%%   Z中位 %.2f   超额中位 %.0f%%\n",
              res[[as.character(dl)]]$detection_rate,
              res[[as.character(dl)]]$median_Z,
              res[[as.character(dl)]]$median_excess_pct))
}

R <- do.call(rbind, res)
write.csv(R, file.path(OUT, "spikein_simulation_corrected.csv"), row.names = FALSE)
cat("\n", strrep("=", 90), "\n Spiked-in 模拟 (修正 P 值) 汇总\n", strrep("=", 90), "\n", sep = "")
print(R, row.names = FALSE)
cat(sprintf("\n  Z 门槛 (BH 秩1, m=%d) = %.2f\n", M_TESTS, z_thresh))
cat("  [saved] spikein_simulation_corrected.csv\n")
