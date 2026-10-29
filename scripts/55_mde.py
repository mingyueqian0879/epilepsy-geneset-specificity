#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 55: 最小可检出富集效应 (MDE) —— 用于审稿自查
  问题: 在 248 个检验 + BH 校正下, 基因集需要多大超额富集才能被检出?
"""
import numpy as np
import pandas as pd

OUT = "out"
R = pd.read_csv(OUT + r"\specificity_matrix.csv")
Rz = R.dropna(subset=["Z"])
S = pd.read_csv(OUT + r"\dataset_summary.csv")

print("=" * 92)
print(" 最小可检出富集效应 (Minimum Detectable Enrichment)")
print("=" * 92)
print("  计算: 给定检验数 m 与 BH 阈值, 求所需 Z; 再换算为基因集超额差异表达基因数\n")
print("  %-12s %-7s %-9s %-10s %-12s %s" % (
    "队列", "检验数", "DEG率", "BH所需Z", "k=100时", "可检出的相对超额"))

for _, row in S.iterrows():
    ds = row.dataset
    sub = Rz[Rz.dataset == ds]
    if not len(sub):
        continue
    m = len(sub)
    # BH: 需 P < 0.05*1/m (最严格的第一个)
    p_need = 0.05 / m
    from scipy import stats
    z_need = stats.norm.isf(p_need / 2)
    p = row.n_deg / row.n_genes
    k = 100
    sd_counts = np.sqrt(k * p * (1 - p))
    exp_counts = k * p
    extra = z_need * sd_counts
    rel = extra / exp_counts if exp_counts > 0 else np.nan
    print("  %-12s %-7d %-9s %-10s %-12s %s" % (
        ds, m, "%.1f%%" % (100 * p), "%.2f" % z_need,
        "expect %.1f" % exp_counts,
        ("+%.0f%% (%.1f extra)" % (100 * rel, extra)) if np.isfinite(rel) else "n.a."))

print()
print("=" * 92)
print(" 解读")
print("=" * 92)
print("  若某个真实基因集只比全基因组背景多 10-15% 的差异表达基因(温和但真实的富集),")
print("  本研究能否检出? 见上表最后一列 —— 所需超额普遍远高于 10-15%。")
print()
print("=" * 92)
print(" 阳性对照的实际表现 (需要多大效应才能被 FDR 检出)")
print("=" * 92)
pos = Rz[Rz.gene_set.str.startswith("POS_")]
for _, r in pos.iterrows():
    need = 0.05 / len(Rz[Rz.dataset == r.dataset])
    from scipy import stats
    z_need = stats.norm.isf(need / 2)
    print("  %-12s %-22s n=%-3d obs=%-3d exp=%-6.2f Z=%+.2f  FDR=%.3f  (需Z>=%.2f才显著)" % (
        r.dataset, r.gene_set, r.n, r.obs, r.exp, r.Z, r.FDR, z_need))
print()
print("  -> 小胶质标志物 Z=+2.79 / +6.50 是最强的信号, 但仍未过 FDR 门槛")
