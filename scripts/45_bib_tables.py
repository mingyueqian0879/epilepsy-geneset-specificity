#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 45: 生成 BIB 投稿包正文表格 (Table 1-3)"""
import os
import numpy as np
import pandas as pd

OUT = "out"
PKG = "05_Tables"
os.makedirs(PKG, exist_ok=True)

S = pd.read_csv(os.path.join(OUT, "dataset_summary.csv"))
LP = pd.read_csv(os.path.join(OUT, "label_permutation_200perm.csv"))
CC = pd.read_csv(os.path.join(OUT, "cell_composition_analysis.csv"))
R = pd.read_csv(os.path.join(OUT, "specificity_matrix.csv"))

INFO = {
 "GSE256068": ("TLE with hippocampal sclerosis vs non-epileptic controls",
               "Temporal cortex / hippocampus", "Illumina NovaSeq 6000"),
 "GSE213488": ("Focal cortical dysplasia IIa/IIb vs autopsy controls",
               "Cortical grey matter", "Illumina HiSeq 2500"),
 "GSE128300": ("Focal cortical dysplasia / HS-epilepsy vs non-epilepsy",
               "Cortex", "Illumina HiSeq 2500"),
 "GSE272142": ("Temporal lobe epilepsy vs controls",
               "Dentate gyrus", "Illumina HiSeq 2500"),
}

# ---------------- Table 1 ----------------
deg = CC[CC.marker.str.startswith("DEG_")].pivot_table(
    index="dataset", columns="marker", values="obs", aggfunc="first")
t1 = S.merge(LP[["dataset", "null_median", "null_p95", "emp_P"]], on="dataset", how="left")
t1["condition"] = t1.dataset.map(lambda d: INFO[d][0])
t1["tissue"] = t1.dataset.map(lambda d: INFO[d][1])
t1["platform"] = t1.dataset.map(lambda d: INFO[d][2])
t1["DEG_modelB"] = t1.dataset.map(deg["DEG_modelB_celladj"])
t1["reduction_pct"] = ((t1.n_deg - t1.DEG_modelB) / t1.n_deg.replace(0, np.nan) * 100).round(1)
t1 = t1[["dataset", "condition", "tissue", "platform", "n_case", "n_ctrl", "n_genes",
         "n_deg", "deg_rate", "DEG_modelB", "reduction_pct",
         "null_median", "null_p95", "emp_P"]]
t1.columns = ["Dataset", "Condition", "Tissue", "Platform", "Cases", "Controls",
              "Genes tested", "DEGs (Model A)", "DEG rate (%)",
              "DEGs (Model B, cell-adjusted)", "Reduction (%)",
              "Null median DEGs", "Null 95th pct", "Empirical P"]
t1.to_csv(os.path.join(PKG, "Table1_datasets_and_calibration.csv"), index=False)
print("[Table 1]")
print(t1.to_string(index=False))

# ---------------- Table 2 ----------------
Rz = R.dropna(subset=["Z"])
rows = []
for ds in ["GSE256068", "GSE213488", "GSE128300", "GSE272142"]:
    sub = Rz[Rz.dataset == ds]
    if not len(sub):
        rows.append({"Cohort": ds, "Gene sets tested": 0, "Median Z": np.nan,
                     "IQR Z": "—", "|Z|>2 (%)": "—", "Nominal P<0.05 (%)": "—",
                     "FDR<0.05 (n)": 0, "Expected at 5% (n)": "—"})
        continue
    rows.append({
        "Cohort": ds, "Gene sets tested": len(sub),
        "Median Z": round(sub.Z.median(), 3),
        "IQR Z": "%.2f to %.2f" % (sub.Z.quantile(.25), sub.Z.quantile(.75)),
        "|Z|>2 (%)": round(100 * (sub.Z.abs() > 2).mean(), 1),
        "Nominal P<0.05 (%)": round(100 * (sub.P < 0.05).mean(), 1),
        "FDR<0.05 (n)": int((sub.FDR < 0.05).sum()),
        "Expected at 5% (n)": round(0.05 * len(sub), 1)})
t2 = pd.DataFrame(rows)
t2.to_csv(os.path.join(PKG, "Table2_geneset_enrichment_summary.csv"), index=False)
print("\n[Table 2]")
print(t2.to_string(index=False))

# ---------------- Table 3 ----------------
hot = Rz[Rz.family == "HOT"].copy()
LABEL = {"Lactylation_LRG332": "Lactylation-related genes (332)",
         "Ferroptosis_core": "Ferroptosis (core)",
         "Cuproptosis_core": "Cuproptosis (core)",
         "Disulfidptosis_core": "Disulfidptosis (core)",
         "Pyroptosis_core": "Pyroptosis (core)",
         "m6A_regulators": "m6A regulators",
         "Autophagy_core": "Autophagy (core)",
         "Glycolysis_core": "Glycolysis (core)",
         "POS_Astrocyte": "Positive control: astrocyte markers",
         "POS_Microglia": "Positive control: microglial markers",
         "POS_Neuron": "Positive control: neuronal markers",
         "POS_Oligodendrocyte": "Positive control: oligodendrocyte markers"}
ORDER = ["Lactylation_LRG332", "Ferroptosis_core", "Cuproptosis_core",
         "Disulfidptosis_core", "Pyroptosis_core", "m6A_regulators",
         "Autophagy_core", "Glycolysis_core",
         "POS_Astrocyte", "POS_Microglia", "POS_Neuron", "POS_Oligodendrocyte"]
rows = []
for g in ORDER:
    sub = hot[hot.gene_set == g]
    if not len(sub):
        continue
    rec = {"Gene set": LABEL.get(g, g), "Genes detected": int(sub.n.max())}
    for ds, tag in [("GSE256068", "GSE256068"), ("GSE213488", "GSE213488"),
                    ("GSE128300", "GSE128300")]:
        s = sub[sub.dataset == ds]
        rec[tag + " Z"] = round(float(s.Z.iloc[0]), 2) if len(s) else "n.t."
        rec[tag + " FDR"] = (round(float(s.FDR.iloc[0]), 3) if len(s) else "n.t.")
    fdr_vals = [float(x) for x in sub.FDR if pd.notna(x)]
    rec["Minimum FDR"] = round(min(fdr_vals), 3) if fdr_vals else np.nan
    rows.append(rec)
t3 = pd.DataFrame(rows)
t3.to_csv(os.path.join(PKG, "Table3_candidate_genesets.csv"), index=False)
print("\n[Table 3]")
print(t3.to_string(index=False))

# ---------------- 补充表复制 ----------------
import shutil
SUP_SRC = "supplementary"
SUP_DST = "06_Supplementary"
os.makedirs(SUP_DST, exist_ok=True)
for f in os.listdir(SUP_SRC):
    if f.endswith((".csv", ".md")):
        shutil.copy2(os.path.join(SUP_SRC, f), os.path.join(SUP_DST, f))
print("\n[补充材料已复制 %d 个文件]" % len(os.listdir(SUP_DST)))

print("\n" + "=" * 70)
print("投稿包表格清单")
print("=" * 70)
for d in [PKG, SUP_DST]:
    print("\n" + d)
    for f in sorted(os.listdir(d)):
        print("   %-46s %7.1f KB" % (f, os.path.getsize(os.path.join(d, f)) / 1024))
