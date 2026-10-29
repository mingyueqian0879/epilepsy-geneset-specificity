# Step 107: 复核手稿的三个 DEG 率
#   (a) TLE-HS(64) vs 全部对照(27)                     手稿称 59.7%
#   (b) TLE-HS(64) vs ControlCortex(14)  无组织协变量   手稿称 44.3%
#   (c) TLE-HS(64) vs ControlCortex(14)  有组织协变量   手稿称  0.5%
suppressMessages({library(limma); library(jsonlite)})
BASE <- Sys.getenv("PROJECT_ROOT", unset = "."); OUT <- file.path(BASE, "out")

expr_df <- read.delim(file.path(OUT, "expr_matrix.tsv"), row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
expr <- as.matrix(expr_df); storage.mode(expr) <- "double"
meta <- read.csv(file.path(OUT, "design_aligned.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$matrix_col
meta <- meta[colnames(expr), ]

run <- function(keep, formula, tag) {
  e <- expr[, keep]; m <- meta[keep, ]
  m$grp <- factor(ifelse(m$disease == "TLE-HS", "Case", "Ctrl"), levels = c("Ctrl", "Case"))
  design <- model.matrix(formula, data = m)
  fit <- eBayes(lmFit(e, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grpCase", number = Inf, sort.by = "P")
  n <- sum(tt$adj.P.Val < 0.05)
  cat(sprintf("%-58s  n=%3d  DEG=%6d / %5d = %5.1f%%\n",
              tag, length(keep), n, nrow(tt), 100 * n / nrow(tt)))
  invisible(tt)
}

cat("===============================================================\n")
cat(" GSE256068 三个模型复核\n")
cat("===============================================================\n")
a <- run(meta$disease %in% c("TLE-HS", "ControlCortex", "ControlHippocampus"),
         ~ grp, "(a) TLE-HS(64) vs 全部对照(27)          无协变量")
cat("\n    注: (a) 若含组织协变量则为:\n")
a2 <- run(meta$disease %in% c("TLE-HS", "ControlCortex", "ControlHippocampus"),
          ~ grp + factor(tissue), "        TLE-HS(64) vs 全部对照(27)          加组织协变量")

cat("\n")
b <- run(meta$disease %in% c("TLE-HS", "ControlCortex"),
         ~ grp, "(b) TLE-HS(64) vs ControlCortex(14)     无组织协变量")
c3 <- run(meta$disease %in% c("TLE-HS", "ControlCortex"),
          ~ grp + factor(tissue), "(c) TLE-HS(64) vs ControlCortex(14)     加组织协变量")

cat("\n")
cat("ControlCortex 的组织构成:\n")
m <- meta[meta$disease == "ControlCortex", ]
print(table(m$tissue))
cat("\nTLE-HS 的组织构成:\n")
print(table(meta$tissue[meta$disease == "TLE-HS"]))
