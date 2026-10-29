# ============================================================
# Step 63: 方法对比 —— 本管线 vs 既有特异性检验方法
#  比较: (a) 表达量匹配置换 (本文)  (b) camera (competitive)
#        (c) roast / fry (self-contained)  (d) 朴素超几何 (领域常用)
#  输出: FDR 显著数、排名相关性、结论一致性
#  注: SetRank 因与 R 4.6.1/Bioc 3.23 不兼容(安装时段错误)未能纳入
# ============================================================
suppressMessages({library(limma); library(edgeR)})
set.seed(20260922)
NPERM <- 2000

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
MOUT <- file.path(BASE, "out", "multi")
OUT  <- file.path(BASE, "out")

## ---------- 基因集 ----------
read_gmt <- function(f) {
  out <- list()
  for (l in readLines(f)) {
    p <- strsplit(l, "\t")[[1]]
    if (length(p) >= 3) out[[p[[1]]]] <- unique(p[-c(1, 2)])
  }
  out
}
MS <- file.path(BASE, "data", "msigdb")
GS <- c(read_gmt(file.path(MS, "h.all.v2024.1.Hs.symbols.gmt")),
        read_gmt(file.path(MS, "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")))
cat("[gene sets]", length(GS), "\n")

mapdf <- read.csv(file.path(MOUT, "ensg2symbol_gencode.csv"), stringsAsFactors = FALSE)
E2S <- setNames(mapdf$symbol, mapdf$ensg)

matched_Z <- function(is_deg, ave, hit, nperm = NPERM) {
  k <- sum(hit); if (k < 5) return(NA_real_)
  obs <- sum(is_deg[hit])
  br <- unique(quantile(ave, probs = seq(0, 1, 0.05), na.rm = TRUE))
  st <- cut(ave, breaks = br, include.lowest = TRUE, labels = FALSE)
  st[is.na(st)] <- 1L
  need <- tabulate(st[hit], nbins = max(st)); ib <- split(seq_along(st), st)
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

compare_one <- function(ds) {
  cat("\n", strrep("=", 84), "\n [", ds, "]\n", sep = "")
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
    design <- model.matrix(~ factor(dd$group, levels = c("Control", "Case")))
    vm <- voom(dge, design, plot = FALSE)
    fit <- eBayes(lmFit(vm, design), trend = TRUE, robust = TRUE)
    EXPR <- vm$E
  } else {
    E <- E[rowSums(is.finite(E)) == ncol(E), , drop = FALSE]
    design <- model.matrix(~ factor(dd$group, levels = c("Control", "Case")))
    fit <- eBayes(lmFit(E, design), trend = TRUE, robust = TRUE)
    EXPR <- E
  }
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  tt$symbol <- E2S[sub("\\..*$", "", rownames(tt))]
  tt <- tt[!is.na(tt$symbol), ]
  is_deg <- tt$adj.P.Val < 0.05
  cat(sprintf("  genes=%d  DEG=%d (%.1f%%)\n", nrow(tt), sum(is_deg), 100 * mean(is_deg)))

  ## 基因集矩阵 (EXPR 行名是 ENSG, tt 行名同)
  syms <- tt$symbol
  idx <- match(rownames(tt), rownames(EXPR))   # 按 ENSG 对齐
  keep2 <- !is.na(idx)
  M <- EXPR[idx[keep2], , drop = FALSE]

  hitlist <- lapply(GS, function(g) which(toupper(syms[keep2]) %in% toupper(g)))
  names(hitlist) <- names(GS)
  ok <- sapply(hitlist, length) >= 5
  hitlist <- hitlist[ok]
  cat(sprintf("  gene sets testable = %d\n", length(hitlist)))

  ## (a) 本文方法
  z_a <- sapply(hitlist, function(h) matched_Z(is_deg[keep2], tt$AveExpr[keep2], 
                                               seq_along(syms[keep2]) %in% h))
  p_a <- 2 * pnorm(-abs(z_a)); fdr_a <- p.adjust(p_a, "BH")

  ## (b) camera (competitive, 校正基因间相关)
  cam <- tryCatch(camera(M, hitlist, design), error = function(e) NULL)
  ## (c) fry (self-contained 快速近似)
  fryr <- tryCatch(fry(M, hitlist, design), error = function(e) NULL)
  ## (d) 朴素超几何
  N <- length(is_deg[keep2]); K <- sum(is_deg[keep2])
  p_d <- sapply(hitlist, function(h) {
    k <- length(h); o <- sum(is_deg[keep2][h])
    phyper(o - 1, K, N - K, k, lower.tail = FALSE)
  })
  fdr_d <- p.adjust(p_d, "BH")

  out <- data.frame(gene_set = names(hitlist), n = sapply(hitlist, length),
                    Z_matched = round(z_a, 2),
                    P_matched = p_a, FDR_matched = fdr_a)
  if (!is.null(cam)) {
    out$P_camera <- cam$PValue[match(out$gene_set, rownames(cam))]
    out$FDR_camera <- cam$FDR[match(out$gene_set, rownames(cam))]
  }
  if (!is.null(fryr)) {
    out$P_fry <- fryr$PValue[match(out$gene_set, rownames(fryr))]
    out$FDR_fry <- fryr$FDR[match(out$gene_set, rownames(fryr))]
  }
  out$P_hyper <- p_d; out$FDR_hyper <- fdr_d

  cat("\n  --- 各方法的 FDR<0.05 显著数 ---\n")
  cat(sprintf("    %-28s %d / %d\n", "本文(表达量匹配置换)", sum(out$FDR_matched < 0.05), nrow(out)))
  if ("FDR_camera" %in% names(out)) cat(sprintf("    %-28s %d / %d\n", "camera (competitive)", sum(out$FDR_camera < 0.05), nrow(out)))
  if ("FDR_fry" %in% names(out))    cat(sprintf("    %-28s %d / %d\n", "fry (self-contained)", sum(out$FDR_fry < 0.05), nrow(out)))
  cat(sprintf("    %-28s %d / %d\n", "朴素超几何(领域常用)", sum(out$FDR_hyper < 0.05), nrow(out)))

  ## 排名相关性
  cat("\n  --- 排名相关性 (Spearman, 基于 P 值) ---\n")
  rk <- function(p) rank(p, na.last = "keep")
  pairs <- list(c("P_matched", "P_camera", "本文 vs camera"),
                c("P_matched", "P_fry", "本文 vs fry"),
                c("P_matched", "P_hyper", "本文 vs 超几何"),
                c("P_camera", "P_fry", "camera vs fry"))
  for (pp in pairs) {
    if (all(pp[1:2] %in% names(out))) {
      ok2 <- is.finite(out[[pp[1]]]) & is.finite(out[[pp[2]]])
      if (sum(ok2) > 5) {
        ct <- cor.test(rk(out[[pp[1]]][ok2]), rk(out[[pp[2]]][ok2]), method = "spearman")
        cat(sprintf("    %-22s rho = %+.3f  (n=%d)\n", pp[3], ct$estimate, sum(ok2)))
      }
    }
  }
  out$dataset <- ds
  out
}

RES <- rbind(compare_one("GSE256068"), compare_one("GSE213488"))
write.csv(RES, file.path(OUT, "method_comparison.csv"), row.names = FALSE)

cat("\n", strrep("=", 84), "\n 跨方法一致性汇总\n", strrep("=", 84), "\n", sep = "")
for (ds in unique(RES$dataset)) {
  s <- RES[RES$dataset == ds, ]
  cat(sprintf("\n[%s] n=%d\n", ds, nrow(s)))
  sets <- intersect(c("FDR_matched","FDR_camera","FDR_fry","FDR_hyper"), names(s))
  lab <- c(FDR_matched="本文", FDR_camera="camera", FDR_fry="fry", FDR_hyper="超几何")
  for (i in seq_along(sets)) for (j in seq_along(sets)) {
    if (i >= j) next
    a <- s[[sets[i]]] < 0.05; b <- s[[sets[j]]] < 0.05
    if (sum(a) + sum(b) == 0) next
    tab <- table(a, b)
    jac <- sum(a & b) / max(sum(a | b), 1)
    cat(sprintf("   %-10s vs %-10s  双方显著交集 %d   Jaccard %.2f\n",
                lab[sets[i]], lab[sets[j]], sum(a & b), jac))
  }
}
cat("\n[saved] method_comparison.csv\n")
