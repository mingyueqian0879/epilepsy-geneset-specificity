# Input data

No data are committed to this repository. Everything is public and is downloaded by the
scripts, but the expected file names are fixed, so the layout below must be respected.

```
data/
├── msigdb/
│   ├── h.all.v2024.1.Hs.symbols.gmt
│   └── c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt
├── GSE256068_processed_data.csv.gz        # log2 TPM, 19,999 genes x 162 samples
├── multi/                                 # count matrices for GSE213488 / GSE128300 / GSE272142
│   ├── GSE213488_*.csv
│   ├── GSE128300_*.csv
│   └── GSE272142_*.csv
├── gencode.v44.annotation.gtf.gz
└── cheng2023_tables/                      # source of the 332 lactylation-related genes
    └── Table S1 List of 332 lactylation-related genes.xlsx
```

## Where each file comes from

| File | Source | How |
|---|---|---|
| `GSE256068_processed_data.csv.gz` | GEO accession GSE256068 | series matrix / supplementary processed data |
| `multi/GSE213488_*.csv`, `multi/GSE128300_*.csv`, `multi/GSE272142_*.csv` | GEO | `python scripts/26_multi_download.py`, `python scripts/27_download_multi.py` produce the raw files; `28`–`30` standardise them into `out/multi/` |
| `msigdb/*.gmt` | MSigDB v2024.1.Hs (Broad Institute) | https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp — download `h.all` and `c2.cp.kegg_legacy`, symbol format. Registration is free; the files are not redistributed here for licensing reasons |
| `gencode.v44.annotation.gtf.gz` | GENCODE v44 | https://www.gencodegenes.org/human/release_44.html |
| `cheng2023_tables/...xlsx` | Cheng et al. 2023, Pharmaceuticals 16:644 (Table S1) | publisher supplementary material |

## Derived inputs

Two files are produced by the scripts themselves and are required downstream:

| File | Produced by | Consumed by |
|---|---|---|
| `out/LRG_final_symbols.txt` | `06`–`09` (lactylation gene list) | `32_core_analysis.R` |
| `out/ensg2symbol_gencode.csv` | `31_gencode_map.py` | `32_core_analysis.R`, `116_ladder_symbol_frame.R` |

## Sanity checks

Before running the pipeline, confirm that `GSE256068_processed_data.csv.gz` decompresses to
19,999 genes × 162 samples (91 cases/controls used for the analysis) and that both GMT
files load as 50 and 186 gene sets respectively:

```bash
Rscript -e 'source("scripts/32_core_analysis.R")'   # stops early if the GMT files are missing
```
