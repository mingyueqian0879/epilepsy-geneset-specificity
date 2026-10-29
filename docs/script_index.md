# Script index

Every script in `scripts/`, with its role and first comment line.
Roles: **core** = produces a number, figure or table in the manuscript;
**data** = acquires/prepares input data; **figures** / **tables** = presentation;
**verification** = re-derives or audits a reported value; **exploratory** =
abandoned or preliminary branches kept for provenance; **manuscript** =
typesetting, reference and word-count utilities (not needed to reproduce results).

## core (20 scripts)

| script | lines | first comment |
|---|---|---|
| `32_core_analysis.R` | 179 | Step 32: 核心分析 —— 多数据集 x 多基因集 特异性检验矩阵 |
| `33_field_survey.py` | 128 | !/usr/bin/env python |
| `34_analyze_matrix.py` | 77 | !/usr/bin/env python |
| `35_label_permutation.R` | 144 | Step 35: 校准检验 —— 标签置换 (label permutation) |
| `35b_labelperm_fast.R` | 78 | Step 35b: 标签置换校准 (快速版) —— 只检验 DEG 数 |
| `36_cell_composition.R` | 119 | Step 36: 细胞组成混杂的多数据集量化 |
| `55_mde.py` | 61 | !/usr/bin/env python |
| `56_p0_numbers.py` | 119 | !/usr/bin/env python |
| `57_fig7_v2.py` | 64 | !/usr/bin/env python |
| `58_suppl_S7.py` | 120 | !/usr/bin/env python |
| `59_p0_verify.py` | 94 | !/usr/bin/env python |
| `60_spikein.R` | 109 | Step 60: Spiked-in 模拟 |
| `61_perm_resolution.py` | 96 | !/usr/bin/env python |
| `62_spikein_fixed.R` | 113 | Step 62: Spiked-in 模拟 (修正版) |
| `63_method_compare.R` | 168 | Step 63: 方法对比 —— 本管线 vs 既有特异性检验方法 |
| `64_camera_check.py` | 56 | !/usr/bin/env python |
| `65_recompute_tables.py` | 90 | !/usr/bin/env python |
| `115_model_ladder.R` | 88 | Step 115: 用脚本 36 的原始标记基因与预处理，重算 GSE256068 的四级模型阶梯 |
| `116_ladder_symbol_frame.R` | 89 | Step 116: 在"基因集分析口径"(有 HGNC symbol, 19958 基因) 上重算完整模型阶梯 |
| `117_rebuild_cellcomp.R` | 97 | Step 117: 在"基因集分析口径"(有 HGNC symbol) 上统一重算 TableS4 与 Table 1 的相关数字 |

## data (15 scripts)

| script | lines | first comment |
|---|---|---|
| `02_screen_datasets.py` | 93 | !/usr/bin/env python |
| `03_enumerate_geo.py` | 90 | !/usr/bin/env python |
| `04_check_targets.py` | 87 | !/usr/bin/env python |
| `05_load_gse256068.py` | 106 | !/usr/bin/env python |
| `06_find_lrg_source.py` | 65 | !/usr/bin/env python |
| `07_trace_lrg.py` | 65 | !/usr/bin/env python |
| `08_build_lrg.py` | 77 | !/usr/bin/env python |
| `09_build_lrg2.py` | 88 | !/usr/bin/env python |
| `25_screen2.py` | 117 | !/usr/bin/env python |
| `26_multi_download.py` | 107 | !/usr/bin/env python |
| `27_download_multi.py` | 61 | !/usr/bin/env python |
| `28_inspect_multi.py` | 40 | !/usr/bin/env python |
| `29_multi_design.py` | 83 | !/usr/bin/env python |
| `30_prep_multi.py` | 204 | !/usr/bin/env python |
| `31_gencode_map.py` | 51 | !/usr/bin/env python |

## figures (18 scripts)

| script | lines | first comment |
|---|---|---|
| `37_figures.py` | 167 | !/usr/bin/env python |
| `38_figures2.py` | 118 | !/usr/bin/env python |
| `39_fig5_redo.py` | 79 | !/usr/bin/env python |
| `40_fig8_suppl.py` | 138 | !/usr/bin/env python |
| `46_graphical_abstract.py` | 102 | !/usr/bin/env python |
| `47_graphabs_v2.py` | 75 | !/usr/bin/env python |
| `66_figs_update.py` | 120 | !/usr/bin/env python |
| `113_sample_level_figure.py` | 122 | !/usr/bin/env python |
| `114_fig2_with_samples.py` | 124 | !/usr/bin/env python |
| `126b_layout_diag.py` | 143 | !/usr/bin/env python |
| `126c_final_params.py` | 129 | !/usr/bin/env python |
| `126d_fig2c_candidates.py` | 150 | !/usr/bin/env python |
| `127_fix_fig2_figS1.py` | 295 | !/usr/bin/env python |
| `127_sync_frontmatter.py` | 111 | !/usr/bin/env python |
| `128_fig2_geometry_check.py` | 143 | !/usr/bin/env python |
| `133_fig8_legend_check.py` | 121 | !/usr/bin/env python |
| `134_fix_fig8.py` | 131 | !/usr/bin/env python |
| `134a_fig8_legend_scan.py` | 120 | !/usr/bin/env python |

## tables (5 scripts)

| script | lines | first comment |
|---|---|---|
| `45_bib_tables.py` | 124 | !/usr/bin/env python |
| `50_bib_tables_fixed.py` | 126 | !/usr/bin/env python |
| `77_suppl_final.py` | 81 | !/usr/bin/env python |
| `119_update_tables.py` | 95 | !/usr/bin/env python |
| `120_wire_S12_rebuild.py` | 51 | !/usr/bin/env python |

## verification (17 scripts)

| script | lines | first comment |
|---|---|---|
| `43_verify_manuscript.py` | 80 | !/usr/bin/env python |
| `79_final_verify.py` | 98 | !/usr/bin/env python |
| `103_final_check_v2.py` | 153 | !/usr/bin/env python |
| `105_sample_audit.py` | 86 | !/usr/bin/env python |
| `106_recheck_deg.py` | 23 | !/usr/bin/env python |
| `107_verify_deg_rates.R` | 47 | Step 107: 复核手稿的三个 DEG 率 |
| `108_hunt_05pct.R` | 79 | Step 108: 穷举搜索能产生 ~0.5% DEG 率的模型 |
| `109_figure_table_crosscheck.py` | 118 | !/usr/bin/env python |
| `110_remaining_checks.py` | 81 | !/usr/bin/env python |
| `111_labelperm_trace.py` | 50 | !/usr/bin/env python |
| `112_survey_power_check.py` | 74 | !/usr/bin/env python |
| `126_fig_legend_diag.py` | 142 | !/usr/bin/env python |
| `126_sync_admin.py` | 98 | !/usr/bin/env python |
| `129_layout_ocr_audit.py` | 134 | !/usr/bin/env python |
| `130_layout_geometry_audit.py` | 223 | !/usr/bin/env python |
| `131_layout_audit_fixed_inputs.py` | 207 | !/usr/bin/env python |
| `132_band_compare.py` | 16 | from PIL import Image |

## exploratory (17 scripts)

| script | lines | first comment |
|---|---|---|
| `01_parse.py` | 102 | !/usr/bin/env python |
| `10_limma.R` | 129 | Step 4: limma 差异表达分析  —— TLE-HS vs Control |
| `11_verify_signal.R` | 98 | Step 5: 关键验证 —— 是真信号还是混杂? |
| `12_diagnostic.R` | 90 | Step 6: 决定性诊断 —— 多基因集富集 + 随机零分布 Z 分数 |
| `13_cell_adjusted.R` | 96 | Step 7: 细胞组成校正后重测 —— 决定项目生死的实验 |
| `14_tle_internal.R` | 115 | Step 8: 挽救性验证 —— TLE-HS 内部: 乳酸化基因 vs 癫痫严重程度 |
| `15_mr_feasibility.py` | 83 | !/usr/bin/env python |
| `16_gtex_check.py` | 50 | !/usr/bin/env python |
| `17_lrg_eqtl_coverage.py` | 92 | !/usr/bin/env python |
| `18_check_gwas.py` | 31 | !/usr/bin/env python |
| `19_mr_instruments.py` | 97 | !/usr/bin/env python |
| `19b_instruments.py` | 59 | !/usr/bin/env python |
| `20_mr_pilot.py` | 143 | !/usr/bin/env python |
| `21_mr_pilot2.py` | 171 | !/usr/bin/env python |
| `22_mr_background.py` | 135 | !/usr/bin/env python |
| `23_pathway_vs_proteome.R` | 89 | Step ②-E: 决定性对照 —— "蛋白质组衍生列表" vs "真实通路基因集" |
| `24_final_machinery_test.R` | 107 | Step ②-F: 最终验证 —— "乳酸化机器"假说的稳健性 |
