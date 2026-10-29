#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 34: 分析特异性矩阵结果"""
import pandas as pd
import numpy as np

OUT = "out"
R = pd.read_csv(OUT + r"\specificity_matrix.csv")
S = pd.read_csv(OUT + r"\dataset_summary.csv")
R = R.dropna(subset=["Z"])

print("=" * 80)
print(" 1. 数据集层面: DEG 率 与 基因集检验结果")
print("=" * 80)
for _, r in S.iterrows():
    sub = R[R.dataset == r.dataset]
    sig = sub[sub.FDR < 0.05]
    print("\n%s : Case=%d Control=%d  基因=%d  DEG=%d (%.1f%%)" % (
        r.dataset, r.n_case, r.n_ctrl, r.n_genes, r.n_deg, r.deg_rate))
    if len(sub):
        print("   基因集检验 n=%d   Z: 中位%+.3f 均值%+.3f SD=%.3f   |Z|>2: %d" % (
            len(sub), sub.Z.median(), sub.Z.mean(), sub.Z.std(), (sub.Z.abs() > 2).sum()))
        print("   FDR<0.05: %d (%.1f%%)   预期假阳性率 5%% -> 期望 %.1f 个" % (
            len(sig), 100 * len(sig) / len(sub), 0.05 * len(sub)))
    else:
        print("   无 DEG -> 无法检验基因集 (Z 全为 NaN)")

print("\n" + "=" * 80)
print(" 2. FDR<0.05 的基因集 (全部)")
print("=" * 80)
sig = R[R.FDR < 0.05].sort_values(["dataset", "FDR"])
if len(sig):
    print(sig[["dataset", "family", "gene_set", "n", "obs", "exp", "Z", "P", "FDR"]]
          .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
else:
    print("  (无)")

print("\n" + "=" * 80)
print(" 3. 各数据集 Z 分数最高的 10 个基因集")
print("=" * 80)
for ds in R.dataset.unique():
    sub = R[R.dataset == ds].nlargest(10, "Z")
    print("\n[%s]" % ds)
    print(sub[["family", "gene_set", "n", "obs", "exp", "Z", "P", "FDR"]]
          .to_string(index=False, float_format=lambda x: f"{x:.4g}"))

print("\n" + "=" * 80)
print(" 4. 阳性对照 vs 疾病热点基因集")
print("=" * 80)
R["cat"] = np.where(R.family == "HOT", "hot", "MSigDB")
pos = R[R.gene_set.str.startswith("POS_")]
hot = R[R.family == "HOT"]
msig = R[R.family != "HOT"]
print("\n阳性对照 (细胞类型标志物):")
print(pos[["dataset", "gene_set", "n", "obs", "exp", "Z", "P", "FDR"]]
      .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
print("\n热点基因集 Z 分数汇总 (按基因集):")
g = hot.groupby("gene_set").agg(n_ok=("Z", "count"), Z_med=("Z", "median"),
                                Z_max=("Z", "max"), minFDR=("FDR", "min")).reset_index()
print(g.sort_values("Z_max", ascending=False).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print("\nMSigDB (Hallmark+KEGG) Z 分数分布:")
print("   全部: n=%d Z中位=%+.3f | Z>2: %d (%.1f%%)  Z<-2: %d (%.1f%%)  FDR<0.05: %d (%.1f%%)" % (
    len(msig), msig.Z.median(), (msig.Z > 2).sum(), 100 * (msig.Z > 2).mean(),
    (msig.Z < -2).sum(), 100 * (msig.Z < -2).mean(),
    (msig.FDR < 0.05).sum(), 100 * (msig.FDR < 0.05).mean()))
print("   双侧 |Z|>2 比例 = %.1f%%  (正态预期 4.6%%)" % (100 * (msig.Z.abs() > 2).mean()))

print("\n" + "=" * 80)
print(" 5. 关键统计: 全领域基因集检验的校准性")
print("=" * 80)
for ds in R.dataset.unique():
    sub = R[R.dataset == ds]
    if not len(sub):
        continue
    obs5 = (sub.P < 0.05).mean() * 100
    print("  %-12s 名义P<0.05 比例 = %5.1f%%  (标称5%%)   |   FDR<0.05 比例 = %4.1f%%" % (
        ds, obs5, (sub.FDR < 0.05).mean() * 100))
