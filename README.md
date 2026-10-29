# Specificity-testing gene-set associations in epilepsy transcriptomics

Analysis code, parameter configurations and environment specifications for the
manuscript *"Specificity-testing gene-set associations in epilepsy transcriptomics: a
calibrated permutation pipeline and reporting standard"* (Briefings in Bioinformatics,
Problem solving protocol).

The pipeline tests whether a candidate gene set carries more differentially expressed
genes than an expression-matched random set of genes — the specificity test that this
work found reported by 9 of 576 surveyed epilepsy bioinformatics papers.

## What this repository contains — and what it does not

| | |
|---|---|
| **Contains** | every analysis script (R and Python), the fixed random seed, permutation counts, thresholds, gene-set definitions, and the small result tables that back the manuscript's figures and supplementary tables |
| **Does not contain** | the expression data. All input data are public and are downloaded by the scripts in `scripts/`; nothing here is patient-identifiable or access-restricted |
| **Not applicable** | this study generated no new sequencing data, so no data-use agreement or controlled-access procedure applies |

## Data sources

| Resource | Accession / version | Used for |
|---|---|---|
| GSE256068 | GEO | temporal lobe epilepsy with hippocampal sclerosis vs non-epileptic controls (64/27) |
| GSE213488 | GEO | focal cortical dysplasia IIa/IIb vs autopsy controls (13/4) |
| GSE128300 | GEO | FCD / HS-epilepsy vs non-epilepsy (20/8) |
| GSE272142 | GEO | temporal lobe epilepsy vs controls (6/6) |
| MSigDB | v2024.1.Hs | 50 Hallmark sets, 186 KEGG legacy pathways |
| GENCODE | v44 | Ensembl → HGNC symbol mapping |
| PubMed | E-utilities, 22 Sep 2026 | the 576-paper field survey |

See [`data/README.md`](data/README.md) for file names, expected locations and download
instructions.

## Requirements

- **R 4.6.1** with `limma` 3.68.4 and `edgeR` 4.10.1
- **Python 3.12** with pandas, numpy, scipy, matplotlib, python-docx, Pillow, pyliftover
- ~8 GB RAM for the expression-matched permutation step

Exact versions: [`environment/requirements.txt`](environment/requirements.txt),
[`environment/R_packages.md`](environment/R_packages.md). Record your own with
`Rscript environment/session_info.R`.

## Quick start

```bash
git clone https://github.com/mingyueqian0879/epilepsy-geneset-specificity.git
cd epilepsy-geneset-specificity

export PROJECT_ROOT="$PWD"          # Windows PowerShell: $env:PROJECT_ROOT = (Get-Location).Path

Rscript environment/session_info.R  # writes environment/sessionInfo.txt
```

Then, in this order (each step writes into `out/`):

```bash
Rscript scripts/32_core_analysis.R      # 248 gene sets x 4 cohorts, 2,000 permutations  (~40 min)
Rscript scripts/35b_labelperm_fast.R    # label permutation, 200 per cohort               (~30 min)
Rscript scripts/36_cell_composition.R   # cell-composition adjustment
python  scripts/33_field_survey.py      # the 576-paper literature survey
Rscript scripts/62_spikein_fixed.R      # spiked-in simulation (7 magnitudes x 30 replicates)
Rscript scripts/63_method_compare.R     # camera / fry / hypergeometric comparison
Rscript scripts/115_model_ladder.R      # four-model adjustment ladder (Table S12)
```

Everything is deterministic given `set.seed(20260922)`.

## How the manuscript maps onto the code

| Manuscript item | Script(s) | Output |
|---|---|---|
| Figure 1; Table S3 | `33_field_survey.py` | `out/field_survey_papers.json` |
| Figure 2A, Table 1, Table S1 | `32_core_analysis.R`, `35b_labelperm_fast.R`, `105_sample_audit.py` | `out/specificity_matrix.csv`, `out/label_permutation_200perm.csv` |
| Figure 2B (patient-level PCA) | `113_sample_level_figure.py`, `114_fig2_with_samples.py` | `out/sample_meta.json` |
| Figure 2C-E; Figures 3, 4 | `34_analyze_matrix.py`, `37_figures.py`, `38_figures2.py`, `39_fig5_redo.py`, `40_fig8_suppl.py` | `figures/` |
| Table 2; Table S2; Table S9 | `32_core_analysis.R`, `34_analyze_matrix.py` | `out/specificity_matrix.csv` |
| Table 3 | `32_core_analysis.R`, `65_recompute_tables.py` | `out/specificity_matrix.csv` |
| Figure 5; Table S10 | `60_spikein.R`, `62_spikein_fixed.R` | `out/spikein_simulation_corrected.csv` |
| Figure 4A-B; Table S4 | `36_cell_composition.R`, `117_rebuild_cellcomp.R` | `out/cell_composition_analysis.csv` |
| Table S11, S11b | `63_method_compare.R`, `64_camera_check.py` | `out/method_comparison.csv` |
| Table S7 (power) | `58_suppl_S7.py` | `out/` |
| Table S8 (MDE) | `55_mde.py` | `out/` |
| Table S12 (model ladder) | `115_model_ladder.R`, `116_ladder_symbol_frame.R` | `out/` |
| six-item reporting standard (Figure 6) | `57_fig7_v2.py` | `figures/` |

Full index of all scripts: [`docs/script_index.md`](docs/script_index.md).
Output-to-manuscript map: [`docs/outputs_map.md`](docs/outputs_map.md).

## Repository layout

```
scripts/        all analysis scripts, numbered in execution order
results/        small result tables that back the manuscript (regenerated by the scripts)
data/           download instructions and expected file names (no data committed)
environment/    R/Python version specifications
docs/           script index, output map, reproducibility notes
```

## Licence and citation

Code released under the MIT licence — see [`LICENSE`](LICENSE). If you use this
pipeline, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)) and, where
relevant, the archived release DOI.

## Contact

Corresponding author: *to be completed in the title page of the manuscript.*
