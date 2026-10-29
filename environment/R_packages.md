# R packages

Tested with **R 4.6.1** (2026).

| Package | Version | Used for |
|---|---|---|
| limma | 3.68.4 | differential expression, camera, fry |
| edgeR | 4.10.1 | TMM normalisation of count data |
| jsonlite | any | JSON I/O in the R helpers |

Install with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"), ask = FALSE)
install.packages("jsonlite")
```

`SetRank` is deliberately not used: under R 4.6.1 / Bioconductor 3.23 it reproducibly
segfaults, which is why its comparison is absent from the manuscript.
