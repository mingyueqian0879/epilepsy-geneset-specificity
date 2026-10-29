#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 65: 用修正后的 P 值重算 Table 2/3 与 Figure 3/4 的数据
  变更: P 由置换零分布估计 Z 后经正态近似得到 (无 1/(n+1) 分辨率限制)
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = "."
OUT = os.path.join(BASE, "out")
PKG = "."
TAB = os.path.join(PKG, "05_Tables")

R = pd.read_csv(os.path.join(OUT, "specificity_matrix.csv"))
Rz = R.dropna(subset=["Z"]).copy()
Rz["P"] = 2 * stats.norm.sf(np.abs(Rz["Z"]))
Rz["FDR"] = np.nan
for ds in Rz.dataset.unique():
    i = Rz.dataset == ds
    Rz.loc[i, "FDR"] = stats.false_discovery_control(Rz.loc[i, "P"].clip(1e-300), method="bh")
Rz.to_csv(os.path.join(OUT, "specificity_matrix_final.csv"), index=False)

print("=" * 92)
print(" 修正 P 值后的结果")
print("=" * 92)
for ds in ["GSE256068", "GSE213488", "GSE128300"]:
    s = Rz[Rz.dataset == ds]
    if not len(s):
        continue
    print("  %-12s m=%-4d  名义P<0.05=%-4d  FDR<0.05=%d" % (
        ds, len(s), int((s.P < 0.05).sum()), int((s.FDR < 0.05).sum())))

# ---------- Table 2 ----------
rows = []
for ds in ["GSE256068", "GSE213488", "GSE128300", "GSE272142"]:
    s = Rz[Rz.dataset == ds]
    if not len(s):
        rows.append({"Cohort": ds, "Gene sets tested": 0, "Median Z": "n.a.",
                     "IQR of Z": "n.a.", "|Z|>2 (%)": "n.a.", "Nominal P<0.05 (%)": "n.a.",
                     "FDR<0.05 (n)": 0, "Expected at 5% (n)": "n.a."})
        continue
    rows.append({"Cohort": ds, "Gene sets tested": len(s),
                 "Median Z": round(float(s.Z.median()), 3),
                 "IQR of Z": "%.2f to %.2f" % (s.Z.quantile(.25), s.Z.quantile(.75)),
                 "|Z|>2 (%)": round(100 * (s.Z.abs() > 2).mean(), 1),
                 "Nominal P<0.05 (%)": round(100 * (s.P < 0.05).mean(), 1),
                 "FDR<0.05 (n)": int((s.FDR < 0.05).sum()),
                 "Expected at 5% (n)": round(0.05 * len(s), 1)})
pd.DataFrame(rows).to_csv(os.path.join(TAB, "Table2_geneset_enrichment_summary.csv"),
                          index=False)
print("\n[Table 2 已更新]")

# ---------- Table 3 ----------
hot = Rz[Rz.family == "HOT"].copy()
LABEL = {"Lactylation_LRG332": "Lactylation-related genes (332)",
         "Ferroptosis_core": "Ferroptosis (core)", "Cuproptosis_core": "Cuproptosis (core)",
         "Disulfidptosis_core": "Disulfidptosis (core)", "Pyroptosis_core": "Pyroptosis (core)",
         "m6A_regulators": "m6A regulators", "Autophagy_core": "Autophagy (core)",
         "Glycolysis_core": "Glycolysis (core)",
         "POS_Astrocyte": "Calibration probe: astrocyte markers",
         "POS_Microglia": "Calibration probe: microglial markers",
         "POS_Neuron": "Calibration probe: neuronal markers",
         "POS_Oligodendrocyte": "Calibration probe: oligodendrocyte markers"}
ORDER = ["Lactylation_LRG332", "Ferroptosis_core", "Cuproptosis_core", "Disulfidptosis_core",
         "Pyroptosis_core", "m6A_regulators", "Autophagy_core", "Glycolysis_core",
         "POS_Microglia", "POS_Oligodendrocyte", "POS_Astrocyte", "POS_Neuron"]
rows = []
for g in ORDER:
    sub = hot[hot.gene_set == g]
    if not len(sub):
        continue
    rec = {"Gene set": LABEL.get(g, g), "Genes detected": int(sub.n.max())}
    for ds, tag in [("GSE256068", "TLE-HS"), ("GSE213488", "FCD")]:
        s = sub[sub.dataset == ds]
        rec[tag + " Z"] = round(float(s.Z.iloc[0]), 2) if len(s) else "n.t."
        rec[tag + " FDR"] = round(float(s.FDR.iloc[0]), 4) if len(s) else "n.t."
    fv = [float(x) for x in sub.FDR if pd.notna(x)]
    rec["Minimum FDR"] = round(min(fv), 4) if fv else "n.a."
    rows.append(rec)
pd.DataFrame(rows).to_csv(os.path.join(TAB, "Table3_candidate_genesets.csv"), index=False)
print("[Table 3 已更新]")
print(pd.DataFrame(rows).to_string(index=False))

# ---------- 供出图的 Z/FDR 数据 ----------
Rz[["dataset", "family", "gene_set", "n", "obs", "exp", "Z", "P", "FDR"]].to_csv(
    os.path.join(OUT, "specificity_final_for_plots.csv"), index=False)
print("\n[saved] specificity_matrix_final.csv / specificity_final_for_plots.csv")
