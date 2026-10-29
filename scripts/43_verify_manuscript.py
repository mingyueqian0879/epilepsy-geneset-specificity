#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 43: 核验稿件中的关键数字"""
import pandas as pd, numpy as np

OUT = "out"
R = pd.read_csv(OUT + r"\specificity_matrix.csv")
S = pd.read_csv(OUT + r"\dataset_summary.csv")
LP = pd.read_csv(OUT + r"\label_permutation_200perm.csv")
CC = pd.read_csv(OUT + r"\cell_composition_analysis.csv")

print("=" * 78)
print(" 稿件数字核验")
print("=" * 78)

tot = S.n_case.sum() + S.n_ctrl.sum()
print("\n[1] 总样本数 = %d  (稿件写 148)" % tot)
for _, r in S.iterrows():
    print("    %-12s %d + %d = %d" % (r.dataset, r.n_case, r.n_ctrl, r.n_case + r.n_ctrl))

print("\n[2] 基因集检验数")
Rz = R.dropna(subset=["Z"])
print("    有 Z 的检验 = %d" % len(Rz))
print("    记录总数(含NaN) = %d  (稿件 TableS2 写 992)" % len(R))
for ds in R.dataset.unique():
    sub = R[R.dataset == ds]; subz = sub.dropna(subset=["Z"])
    if not len(subz):
        print("    %-12s Z检验=0  总记录=%d" % (ds, len(sub))); continue
    n = len(subz)
    print("    %-12s Z检验=%-4d  名义P<0.05=%5.1f%%  FDR<0.05=%d (%.1f%%)" % (
        ds, n, 100 * (subz.P < 0.05).mean(), (subz.FDR < 0.05).sum(), 100 * (subz.FDR < 0.05).mean()))

print("\n[3] 热点基因集最小 FDR")
hot = Rz[Rz.family == "HOT"]
g = hot.groupby("gene_set").agg(n=("Z", "count"), Zmin=("Z", "min"), Zmax=("Z", "max"),
                                minFDR=("FDR", "min")).reset_index()
print(g.sort_values("minFDR").to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print("    -> 全体最小 FDR = %.3f" % g.minFDR.min())

print("\n[4] 阳性对照 Z")
pos = Rz[Rz.gene_set.str.startswith("POS_")]
print(pos[["dataset", "gene_set", "n", "obs", "exp", "Z", "P"]].to_string(index=False, float_format=lambda x: f"{x:.3g}"))

print("\n[5] 标签置换")
print(LP[["dataset", "obs_deg", "obs_pct", "null_median", "null_p95", "emp_P"]].to_string(index=False))

print("\n[6] 细胞组成")
deg = CC[CC.marker.str.startswith("DEG_")]
p = deg.pivot_table(index="dataset", columns="marker", values="obs", aggfunc="first")
p["降幅%"] = (100 * (p.DEG_modelA - p.DEG_modelB_celladj) / p.DEG_modelA.replace(0, np.nan)).round(1)
print(p.to_string())

print("\n[7] GSE213488 中 FDR<0.05 的全部基因集")
sig = Rz[(Rz.dataset == "GSE213488") & (Rz.FDR < 0.05)].sort_values("FDR")
print(sig[["gene_set", "n", "obs", "exp", "Z", "P", "FDR"]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
print("    炎症/缺氧相关 = %d / %d" % (
    sig.gene_set.str.contains("INFLAMMAT|IL6|IL2|HYPOXIA|CYTOKINE|ANTIGEN", case=False).sum(), len(sig)))

print("\n[8] GSE256068 top Z 基因集")
t = Rz[Rz.dataset == "GSE256068"].nlargest(8, "Z")
print(t[["gene_set", "n", "obs", "exp", "Z", "P", "FDR"]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

print("\n[9] 跨队列相关性")
wide = Rz.pivot_table(index="gene_set", columns="dataset", values="Z", aggfunc="first")
pairs = [("GSE256068", "GSE213488"), ("GSE256068", "GSE128300"), ("GSE213488", "GSE128300")]
rs = []
for a, b in pairs:
    if a in wide and b in wide:
        sub = wide[[a, b]].dropna()
        if len(sub) >= 5:
            r = np.corrcoef(sub[a], sub[b])[0, 1]; rs.append(r)
            print("    %s vs %s : r=%+.3f (n=%d)" % (a, b, r, len(sub)))
print("    mean r = %+.3f  (稿件写 -0.01)" % np.mean(rs))

print("\n[10] 领域调查")
import json
SV = json.load(open(OUT + r"\field_survey_papers.json", encoding="utf-8"))
n = len(SV); nsp = sum(1 for r in SV if r["has_specificity_test"])
print("    论文数=%d  有特异性检验=%d (%.1f%%)  无=%d (%.1f%%)" % (
    n, nsp, 100 * nsp / n, n - nsp, 100 * (n - nsp) / n))
