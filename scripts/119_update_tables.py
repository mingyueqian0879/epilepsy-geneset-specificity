#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 119: 表格更新
 P1-2  TableS6 增列两次运行的 null 95th percentile
 P0-2  Table 1 / TableS1 增列脑区匹配度
 新增  TableS12 四级模型阶梯（正文 §3.9 的数字来源）
"""
import io, os
import pandas as pd

PKG = "."
SUP = os.path.join(PKG, "06_Supplementary")
TAB = os.path.join(PKG, "05_Tables")
NEW = os.path.join(PKG, "06_Supplementary")

# ---------- 1) TableS6：补 null p95 ----------
old = pd.read_csv(os.path.join(SUP, "TableS6_label_permutation_robustness.csv"))
p100 = pd.read_csv("out/label_permutation.csv")
p200 = pd.read_csv("out/label_permutation_200perm.csv")
m100 = dict(zip(p100.dataset, p100.null_deg_p95))
m200 = dict(zip(p200.dataset, p200.null_p95))
rows = []
for _, r in old.iterrows():
    c = r["Cohort"]
    rows.append({
        "Cohort": c,
        "DEGs (100 perm)": int(r["DEGs (100 perm run)"]),
        "Null 95th pct (100 perm)": round(float(m100.get(c, float("nan"))), 2),
        "Empirical P (100 perm)": round(float(r["Empirical P (100 perm)"]), 6),
        "DEGs (200 perm)": int(r["DEGs (200 perm run)"]),
        "Null 95th pct (200 perm)": round(float(m200.get(c, float("nan"))), 2),
        "Empirical P (200 perm)": round(float(r["Empirical P (200 perm)"]), 6)})
S6 = pd.DataFrame(rows)
S6.to_csv(os.path.join(SUP, "TableS6_label_permutation_robustness.csv"), index=False)
print("TableS6 已更新:")
print(S6.to_string(index=False))

# ---------- 2) Table 1 / TableS1：加脑区匹配度 ----------
MATCH = {"GSE256068": "4 / 27", "GSE213488": "4 / 4", "GSE128300": "8 / 8", "GSE272142": "6 / 6"}
NOTE = {"GSE256068": "cases all temporal lobe; controls 13 hippocampus, 8 frontal, 2 parietal, 4 temporal",
        "GSE213488": "all cortical grey matter",
        "GSE128300": "all cortex",
        "GSE272142": "all dentate gyrus"}

T1 = pd.read_csv(os.path.join(TAB, "Table1_datasets_and_calibration.csv"))
T1.insert(3, "Controls region-matched to cases", T1["Dataset"].map(MATCH))
T1.insert(4, "Region note", T1["Dataset"].map(NOTE))
T1.to_csv(os.path.join(TAB, "Table1_datasets_and_calibration.csv"), index=False)
print("\nTable 1 已更新（新增 2 列）:")
print(T1.iloc[:, :5].to_string(index=False))

S1 = pd.read_csv(os.path.join(SUP, "TableS1_datasets.csv"))
S1.insert(4, "controls_region_matched", S1["dataset"].map(MATCH))
S1.insert(5, "region_note", S1["dataset"].map(NOTE))
S1.to_csv(os.path.join(SUP, "TableS1_datasets.csv"), index=False)
print("\nTableS1 已更新（新增 2 列）")
print(S1.iloc[:, :6].to_string(index=False))

# ---------- 3) TableS12：模型阶梯 ----------
LAD = pd.DataFrame([
    {"Comparison": "GSE256068: 64 cases vs all 27 controls", "Model": "M1 group only",
     "Covariates": "none", "DEGs": 11923, "Genes tested": 19999, "DEG rate (%)": 59.62},
    {"Comparison": "GSE256068: 64 cases vs all 27 controls", "Model": "M2 group + tissue",
     "Covariates": "tissue source", "DEGs": 8575, "Genes tested": 19999, "DEG rate (%)": 42.88},
    {"Comparison": "GSE256068: 64 cases vs all 27 controls", "Model": "M3 group + cell scores",
     "Covariates": "6 cell-type scores", "DEGs": 4367, "Genes tested": 19999, "DEG rate (%)": 21.84},
    {"Comparison": "GSE256068: 64 cases vs all 27 controls", "Model": "M4 group + tissue + cell scores",
     "Covariates": "tissue source + 6 cell-type scores", "DEGs": 79, "Genes tested": 19999,
     "DEG rate (%)": 0.40},
    {"Comparison": "GSE256068: 64 cases vs 14 non-epileptic cortex controls", "Model": "M1 group only",
     "Covariates": "none", "DEGs": 12268, "Genes tested": 19999, "DEG rate (%)": 61.34},
    {"Comparison": "GSE256068: 64 cases vs 14 non-epileptic cortex controls", "Model": "M2 group + tissue",
     "Covariates": "tissue source", "DEGs": 8866, "Genes tested": 19999, "DEG rate (%)": 44.33},
    {"Comparison": "GSE256068: 64 cases vs 14 non-epileptic cortex controls", "Model": "M3 group + cell scores",
     "Covariates": "6 cell-type scores", "DEGs": 1907, "Genes tested": 19999, "DEG rate (%)": 9.54},
    {"Comparison": "GSE256068: 64 cases vs 14 non-epileptic cortex controls",
     "Model": "M4 group + tissue + cell scores", "Covariates": "tissue source + 6 cell-type scores",
     "DEGs": 192, "Genes tested": 19999, "DEG rate (%)": 0.96},
])
LAD.to_csv(os.path.join(NEW, "TableS12_model_ladder.csv"), index=False)
print("\n新增 TableS12_model_ladder.csv:")
print(LAD.to_string(index=False))

# ---------- 4) xlsx 同步 ----------
try:
    for d, f in [(TAB, "Table1_datasets_and_calibration"),
                 (SUP, "TableS1_datasets"),
                 (SUP, "TableS6_label_permutation_robustness"),
                 (SUP, "TableS12_model_ladder")]:
        df = pd.read_csv(os.path.join(d, f + ".csv"))
        with pd.ExcelWriter(os.path.join(d, f + ".xlsx"), engine="openpyxl") as w:
            df.to_excel(w, index=False, sheet_name="Data")
        print("  xlsx 同步:", f)
except Exception as e:
    print("  xlsx 同步失败:", str(e)[:90])
