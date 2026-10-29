# ============================================================
# Step 60: Spiked-in 模拟
#  目的: 验证 (a) 零分布在校准状态下的假阳性率; (b) MDE 预测是否准确
#  做法: 在真实 GSE256068 数据中, 向随机选取的 100 基因集注入已知幅度的表达位移,
#        重跑完整流程 (limma + 表达量匹配置换), 记录检出率
# ============================================================
suppressMessages({library(limma); library(edgeR)})
set.seed(20260922)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
MOUT <- file.path(BASE, "out", "multi")
OUT  <- file.path(BASE, "out")

## ---------- 数据 ----------
E <- as.matrix(read.delim(file.path(MOUT, "GSE256068_expr.tsv"),
                          row.names = 1, check.names = FALSE))
dd <- read.csv(file.path(MOUT, "GSE256068_design.csv"), encoding = "UTF-8",
               stringsAsFactors = FALSE)
dd <- dd[dd$sample %in% colnames(E), ]
E <- E[, dd$sample, drop = FALSE]
grp <- factor(dd$group, levels = c("Control", "Case"))
cat(sprintf("[data] %d genes x %d samples  (Case=%d, Control=%d)\n",
            nrow(E), ncol(E), sum(grp == "Case"), sum(grp == "Control")))

SETSIZE  <- 100
NREP     <- 30
NPERM    <- 2000
M_TESTS  <- 248                      # 与正文一致, 用于 BH 阈值
P_THRESH <- 0.05 / M_TESTS           # 秩 1 情形的 BH 门槛

matched_test <- function(is_deg, ave, hit, nperm = NPERM) {
  k <- sum(hit)
  if (k < 5) return(c(obs = NA, exp = NA, Z = NA, P = NA))
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
  c(obs = obs, exp = mean(sims), Z = (obs - mean(sims)) / sd(sims),
    P = (sum(sims >= obs) + 1) / (nperm + 1))
}

run_pipeline <- function(Emat) {
  design <- model.matrix(~ grp)
  fit <- eBayes(lmFit(Emat, design), trend = TRUE, robust = TRUE)
  topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
}

## ---------- 基线 (无注入) ----------
cat("\n[baseline] 无注入, 估计假阳性率\n")
tt0 <- run_pipeline(E)
is0 <- tt0$adj.P.Val < 0.05
cat(sprintf("  genome-wide DEG = %d (%.1f%%)\n", sum(is0), 100 * mean(is0)))

DELTAS <- c(0, 0.15, 0.25, 0.35, 0.5, 0.75)
res <- list()

for (dl in DELTAS) {
  cat(sprintf("\n[spike delta = %.2f] %d 次重复 ...\n", dl, NREP))
  det <- 0; zs <- c(); ps <- c(); excess <- c(); obs_deg <- c()
  for (r in seq_len(NREP)) {
    set.seed(20260922 + r)
    spiked <- sample(rownames(E), SETSIZE)
    Es <- E
    if (dl > 0) {
      case_idx <- which(grp == "Case")
      Es[spiked, case_idx] <- Es[spiked, case_idx] + dl
    }
    tt <- run_pipeline(Es)
    isd <- tt$adj.P.Val < 0.05
    hit <- rownames(tt) %in% spiked
    mt <- matched_test(isd, tt$AveExpr, hit)
    if (!is.na(mt["P"])) {
      ps <- c(ps, mt["P"]); zs <- c(zs, mt["Z"])
      obs_deg <- c(obs_deg, mt["obs"])
      excess <- c(excess, (mt["obs"] - mt["exp"]) / max(mt["exp"], 0.5))
      if (mt["P"] < P_THRESH) det <- det + 1
    }
    if (r %% 10 == 0) cat(sprintf("    rep %d/%d  累计检出 %d\n", r, NREP, det))
  }
  res[[as.character(dl)]] <- data.frame(
    delta = dl, n = length(ps),
    detection_rate = round(100 * det / length(ps), 1),
    median_Z = round(median(zs), 2),
    median_P = signif(median(ps), 3),
    median_excess_pct = round(100 * median(excess), 1),
    median_obs_deg = round(median(obs_deg), 1),
    fp_rate_at_alpha = round(100 * mean(ps < 0.05), 1))
  cat(sprintf("  -> 检出率 %.1f%%  Z中位 %.2f  实测超额中位 %.0f%%\n",
              res[[as.character(dl)]]$detection_rate,
              res[[as.character(dl)]]$median_Z,
              res[[as.character(dl)]]$median_excess_pct))
}

R <- do.call(rbind, res)
write.csv(R, file.path(OUT, "spikein_simulation.csv"), row.names = FALSE)
cat("\n", strrep("=", 88), "\n Spiked-in 模拟汇总\n", strrep("=", 88), "\n", sep = "")
print(R, row.names = FALSE)
cat(sprintf("\n  BH 门槛 (秩1, m=%d) = P < %.2e\n", M_TESTS, P_THRESH))
cat("  [saved] spikein_simulation.csv\n")
