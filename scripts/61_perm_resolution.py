#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 61: 诊断置换分辨率问题
  经验 P 的下限 = 1/(n_perm+1) = 0.0005 (2000 次)
  而 BH 门槛 (秩1, m=248) = 0.05/248 = 2.02e-4
  -> 单个基因集无论多强都无法达到 FDR<0.05
  检验: 改用 Z 的正态近似 P 值后, 结论是否变化?
"""
import numpy as np
import pandas as pd
from scipy import stats

OUT = "out"
R = pd.read_csv(OUT + r"\specificity_matrix.csv")
Rz = R.dropna(subset=["Z"]).copy()

print("=" * 92)
print(" 1) 置换分辨率 vs BH 门槛")
print("=" * 92)
NPERM = 2000
print("  经验 P 下限            = 1/(%d+1) = %.5f" % (NPERM, 1 / (NPERM + 1)))
for m in (248, 248, 96):
    print("  BH 秩1 门槛 (m=%-4d)   = 0.05/%d = %.2e" % (m, m, 0.05 / m))
print()
print("  -> 单个基因集的最大可能经验 P (0.0005) 大于 BH 门槛 (2.02e-4)")
print("  -> 只有当多个基因集并列在 P 下限时, BH 才可能给出 FDR<0.05")
print("     例: 6 个集并列 P=0.0005 -> FDR = 0.0005*248/6 = 0.0207 < 0.05")

print()
print("=" * 92)
print(" 2) 改用 Z 的正态近似 P 值 (无分辨率限制)")
print("=" * 92)
Rz["P_emp"] = Rz["P"]
Rz["P_norm"] = 2 * stats.norm.sf(np.abs(Rz["Z"]))
for ds in Rz.dataset.unique():
    sub = Rz[Rz.dataset == ds].copy()
    m = len(sub)
    for col, lab in [("P_emp", "经验P"), ("P_norm", "正态近似P")]:
        fdr = stats.false_discovery_control(sub[col].clip(1e-300), method="bh")
        n_sig = int((fdr < 0.05).sum())
        print("  %-12s m=%-4d %-10s  名义P<0.05: %-4d   FDR<0.05: %d" % (
            ds, m, lab, int((sub[col] < 0.05).sum()), n_sig))
    print()

print("=" * 92)
print(" 3) 正态近似下各队列 FDR<0.05 的基因集")
print("=" * 92)
for ds in Rz.dataset.unique():
    sub = Rz[Rz.dataset == ds].copy()
    if not len(sub):
        continue
    sub["P_norm"] = 2 * stats.norm.sf(np.abs(sub["Z"]))
    sub["FDR_norm"] = stats.false_discovery_control(sub["P_norm"].clip(1e-300), method="bh")
    sig = sub[sub.FDR_norm < 0.05].sort_values("FDR_norm")
    print("\n  [%s]  m=%d  显著 %d 个" % (ds, len(sub), len(sig)))
    if len(sig):
        for _, r in sig.head(12).iterrows():
            print("      %-46s Z=%+.2f  P=%.2e  FDR=%.4f" % (
                r.gene_set[:46], r.Z, r.P_norm, r.FDR_norm))
    else:
        print("      (无)")

print()
print("=" * 92)
print(" 4) 热点基因集在正态近似下的表现")
print("=" * 92)
hot = Rz[Rz.family == "HOT"].copy()
hot["P_norm"] = 2 * stats.norm.sf(np.abs(hot["Z"]))
rows = []
for ds in hot.dataset.unique():
    sub = hot[hot.dataset == ds].copy()
    sub["FDR_norm"] = stats.false_discovery_control(sub["P_norm"].clip(1e-300), method="bh")
    for _, r in sub.iterrows():
        if r.gene_set.startswith("POS_"):
            continue
        rows.append({"dataset": ds, "gene_set": r.gene_set, "n": int(r.n),
                     "obs": int(r.obs), "exp": round(r.exp, 2),
                     "Z": round(r.Z, 2), "P_norm": r.P_norm,
                     "FDR_norm": round(r.FDR_norm, 4)})
h = pd.DataFrame(rows)
print(h.to_string(index=False, float_format=lambda x: "%.4g" % x))

# 阳性对照
print()
print("  阳性对照 (正态近似 FDR):")
pos = hot[hot.gene_set.str.startswith("POS_")].copy()
for ds in pos.dataset.unique():
    sub = pos[pos.dataset == ds].copy()
    sub["FDR_norm"] = stats.false_discovery_control(
        2 * stats.norm.sf(np.abs(sub["Z"])), method="bh")
    for _, r in sub.sort_values("Z", ascending=False).iterrows():
        print("      %-12s %-22s Z=%+.2f  FDR_norm=%.4f" % (ds, r.gene_set, r.Z, r.FDR_norm))

Rz.to_csv(OUT + r"\specificity_matrix_withPnorm.csv", index=False)
print("\n[saved] specificity_matrix_withPnorm.csv")
