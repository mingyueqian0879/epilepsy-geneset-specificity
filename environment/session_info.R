# Record the exact environment used for a run.
# Usage: Rscript environment/session_info.R
out <- file.path("environment", "sessionInfo.txt")
sink(out)
cat("R version:", R.version.string, "\n\n")
cat("Platform:", R.version$platform, "\n\n")
for (p in c("limma", "edgeR", "jsonlite")) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
  cat(sprintf("%-10s %s\n", p, v))
}
cat("\n")
print(sessionInfo())
sink()
cat("wrote", out, "\n")
