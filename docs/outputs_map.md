# Output → manuscript map

Every number, figure and table in the manuscript traces to one file under `out/`.

| Output file | Produced by | Used in |
|---|---|---|
| `out/specificity_matrix.csv` | `32_core_analysis.R` | Table 2, Table 3, Table S2, Figures 2C–E, 3, 4 |
| `out/specificity_matrix_withPnorm.csv` | `61_perm_resolution.py` | parametric *P* values (replaces the empirical tail) |
| `out/label_permutation_200perm.csv` | `35b_labelperm_fast.R` | Figure 2D–E, Table S5 |
| `out/label_permutation_100perm.csv` | `35b_labelperm_fast.R` (100 permutations) | Table S6 (robustness) |
| `out/cell_composition_analysis.csv` | `36_cell_composition.R` | Figure 4A–B, Table S4 |
| `out/cell_composition_analysis_symbolframe.csv` | `117_rebuild_cellcomp.R` | Table S4 (symbol frame) |
| `out/spikein_simulation_corrected.csv` | `62_spikein_fixed.R` | Figure 5, Table S10 |
| `out/method_comparison.csv` | `63_method_compare.R` | Table S11, S11b |
| `out/field_survey_papers.json` | `33_field_survey.py` | Figure 1, Table S3 |
| `out/sample_meta.json` | `113_sample_level_figure.py` | Figure 2B, Figure S1 |
| `out/multi/` | `28`–`30`, `117` | model ladder (Table S12) |

Machine-readable versions of Tables S1–S12 are supplied with the manuscript as CSV and
XLSX; the copies in `results/` are the same files and exist so that the numbers can be
checked without re-running the pipeline.
