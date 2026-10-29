#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 110: 剩余图表 (Fig5/6/8/9) 与数值一致性细查"""
import io, json, os, re
import numpy as np
import pandas as pd

BASE = "."
OUT = os.path.join(BASE, "out")
PKG = "."
md = io.open(os.path.join(PKG, "03_Manuscript", "Manuscript_BIBformat.md"), encoding="utf-8").read()
body = md[:md.find("## References")]

print("=" * 96)
print(" 1. GSE256068 DEG 计数：三个来源对照")
print("=" * 96)
S = pd.read_csv(os.path.join(OUT, "dataset_summary.csv"))
T1 = pd.read_csv(os.path.join(PKG, "05_Tables", "Table1_datasets_and_calibration.csv"))
CC = pd.read_csv(os.path.join(OUT, "cell_composition_analysis.csv"))
for nm, n, d in [("dataset_summary.csv", int(S[S.dataset == "GSE256068"].n_genes.iloc[0]),
                  int(S[S.dataset == "GSE256068"].n_deg.iloc[0])),
                 ("Table1 (投稿包)", int(T1[T1.Dataset == "GSE256068"]["Genes tested (gene-set analysis)"].iloc[0]),
                  int(T1[T1.Dataset == "GSE256068"]["DEGs, Model A"].iloc[0])),
                 ("cell_composition_analysis.csv (ModelA)",
                  int(CC[(CC.dataset == "GSE256068") & (CC.marker == "DEG_modelA")].n.iloc[0]),
                  int(CC[(CC.dataset == "GSE256068") & (CC.marker == "DEG_modelA")].obs.iloc[0]))]:
    print("  %-42s genes=%6d  DEG=%6d  rate=%.2f%%" % (nm, n, d, 100 * d / n))
print("  手稿正文:                                    genes= 19958  DEG= 11923  rate=%.2f%%" % (100 * 11923 / 19958))
print("  -> 11,923 出自 19,999 基因矩阵；19,958 是基因集分析的过滤后基因数。")
print("     两者混用：正确配对应为 11,905/19,958=59.65% 或 11,923/19,999=59.62%。")

print("\n" + "=" * 96)
print(" 2. GSE213488 DEG 计数 262 vs 263")
print("=" * 96)
print("  dataset_summary n_deg =", int(S[S.dataset == "GSE213488"].n_deg.iloc[0]))
LP = pd.read_csv(os.path.join(OUT, "label_permutation_200perm.csv"))
print("  label_permutation.csv:")
print(LP.to_string(index=False))
print("  手稿正文同时出现:", re.findall(r"26[23]", body)[:12])

print("\n" + "=" * 96)
print(" 3. Fig5 label permutation 数值核对")
print("=" * 96)
for _, r in LP.iterrows():
    print("  %-11s obs_pct=%.2f  null_med=%s null_p95=%s  empP=%s" %
          (r.get("dataset", "?"), r.get("obs_deg_pct", float("nan")),
           r.get("null_median_deg", r.get("null_median", "?")),
           r.get("null_deg_p95", "?"), r.get("emp_P", r.get("empP", "?"))))

print("\n" + "=" * 96)
print(" 4. Fig8 相关性核对")
print("=" * 96)
R = pd.read_csv(os.path.join(OUT, "specificity_matrix_final.csv"))
piv = R.dropna(subset=["Z"]).pivot_table(index="gene_set", columns="dataset", values="Z")
print("  可比较的数据集:", list(piv.columns))
if "GSE256068" in piv and "GSE213488" in piv:
    both = piv[["GSE256068", "GSE213488"]].dropna()
    r = np.corrcoef(both.GSE256068, both.GSE213488)[0, 1]
    print("  GSE256068 vs GSE213488  r = %+.3f  (n=%d)" % (r, len(both)))
    print("  手稿称 r = +0.17 ->", "一致" if abs(r - 0.17) < 0.02 else "** 不一致 **")

print("\n" + "=" * 96)
print(" 5. Fig9 spike-in 数值")
print("=" * 96)
SP = pd.read_csv(os.path.join(OUT, "spikein_simulation_corrected.csv"))
print(SP.to_string(index=False))
print("\n  手稿称: δ=0 假阳性率 0%，+36% 检出率 90%，中位 Z=+0.18")
z0 = SP.iloc[0]
print("  源数据 δ=0 行: 检出率=%s  中位Z=%s" % (z0.get("detection_rate"), z0.get("median_Z")))

print("\n" + "=" * 96)
print(" 6. MDE 表")
print("=" * 96)
M = pd.read_csv(os.path.join(PKG, "06_Supplementary", "TableS8_minimum_detectable_enrichment.csv"))
print(M.to_string(index=False))

print("\n" + "=" * 96)
print(" 7. 手稿中所有 'X% to Y%' 形式的降幅表述")
print("=" * 96)
for m in re.finditer(r"[^.]*\b\d+(?:\.\d+)?%\s*(?:to|\u2192)\s*\d+(?:\.\d+)?%[^.]*\.", body):
    print("  -", m.group().strip()[:170])
