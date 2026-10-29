# Reproducibility notes

## Determinism

Every randomised step uses the seed `20260922`: the expression-matched permutation test
(2,000 replicates per gene set), the label permutation (200 replicates per cohort) and the
spiked-in simulation (30 replicates per spike magnitude, seven magnitudes).

## Two gene universes (do not interchange them)

| Universe | Genes, GSE256068 | Genes, GSE213488 | Used for |
|---|---|---|---|
| full filtered set | 19,999 | 17,139 | label permutation, cell-composition and tissue-adjustment models (Table S12) |
| gene-set analysis set (HGNC symbol available) | 19,958 | 17,010 | all gene-set enrichment tests; DEG rates in Table 1 and Figure 2 |

DEG counts therefore differ slightly between the two frames (11,923 of 19,999 vs 11,905 of
19,958 in GSE256068). Each count is always paired with its own denominator.

## Known trap: `out/label_permutation.csv`

`out/label_permutation.csv` is written by `35_label_permutation.R` and, if that script is
re-run, may hold the **100-permutation** results, whose columns are
`null_deg_median / null_deg_p95 / p_deg`. The manuscript reports the **200-permutation**
results, which live in `out/label_permutation_200perm.csv` (columns
`null_median / null_p95 / emp_P`) and are written by `35b_labelperm_fast.R`.

Every script that needs the 200-permutation values reads that file explicitly —
`38_figures2.py`, `39_fig5_redo.py`, `40_fig8_suppl.py`, `43_verify_manuscript.py`,
`45_bib_tables.py` and `110_remaining_checks.py`. If you regenerate the 100-permutation
file, nothing downstream will silently pick it up.

## Empirical vs parametric *P* values

With 2,000 permutations the smallest attainable empirical *P* is 5 × 10⁻⁴, while the
rank-1 Benjamini–Hochberg threshold for 248 tests is 2.0 × 10⁻⁴: a strong gene set could
never be declared significant from the empirical tail. The pipeline therefore uses the
permutation null only for its mean and standard deviation and takes the tail probability
from the normal approximation `P = 2Φ(-|Z|)` (`61_perm_resolution.py`), a choice validated
by the spiked-in simulation.

## Runtime

| Step | Approximate runtime (8-core desktop) |
|---|---|
| `32_core_analysis.R` | 40 min |
| `35b_labelperm_fast.R` | 30 min |
| `62_spikein_fixed.R` | 45 min |
| `63_method_compare.R` | 10 min |
| everything else | < 5 min |

## Verifying a reported number

`43_verify_manuscript.py`, `59_p0_verify.py` and `79_final_verify.py` re-derive the
manuscript's headline numbers from `out/` and print them for comparison:
148 samples, 592 valid Z scores, 0/248 and 21/248 significant sets, smallest candidate-set
FDR 0.24, label-permutation *P* 0.005 / 0.030 / 0.045 / 1.000, 63% DEG reduction after
cell adjustment, *r* = +0.17 between cohorts, 9/576 papers with a specificity test,
MDE +31% / +297%, 0% false positives and 90% detection at +36% excess, and 0/5/22/148
significant sets for the four specificity tests.
