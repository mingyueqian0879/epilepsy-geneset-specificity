#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 105: 样本水平审计 —— 图 2/图 3 未展示的底层数据到底长什么样
检查：① 分组与组织是否混杂 ② 是否有离群样本 ③ PCA 上组间是否真分离
④ 59.7% 的 DEG 率是否由少数样本主导
"""
import os, pickle
import numpy as np
import pandas as pd

BASE = "."
M = os.path.join(BASE, "out", "multi")

for ds in ["GSE256068", "GSE213488"]:
    X = pd.read_pickle(os.path.join(M, "%s_expr.pkl" % ds))
    D = pd.read_csv(os.path.join(M, "%s_design.csv" % ds))
    print("=" * 88)
    print(" %s   表达矩阵 %s" % (ds, X.shape))
    print("=" * 88)
    # 对齐
    common = [s for s in D["sample"] if s in X.columns]
    D = D[D["sample"].isin(common)].reset_index(drop=True)
    X = X[common]
    print("  样本对齐后: %d" % len(common))

    print("\n  [1] 分组 × 组织 交叉表（混杂检查）")
    ct = pd.crosstab(D["detail"], D["tissue"])
    print(ct.to_string().replace("\n", "\n  "))

    # ---- PCA ----
    print("\n  [2] PCA（前 5 主成分方差占比）")
    Mx = X.values.astype(float)
    Mx = Mx - Mx.mean(axis=1, keepdims=True)
    sd = Mx.std(axis=1, keepdims=True); sd[sd == 0] = 1
    Mx = Mx / sd
    Mx = np.nan_to_num(Mx)
    U, S, Vt = np.linalg.svd(Mx, full_matrices=False)
    var = (S ** 2) / (S ** 2).sum()
    print("    " + "  ".join("PC%d=%.1f%%" % (i + 1, v * 100) for i, v in enumerate(var[:5])))
    pcs = Vt[:5].T * S[:5]          # 样本坐标
    pc = pd.DataFrame(pcs, columns=["PC%d" % (i + 1) for i in range(5)])
    pc["group"] = D["group"].values
    pc["detail"] = D["detail"].values
    pc["tissue"] = D["tissue"].values
    pc["sample"] = D["sample"].values

    print("\n  [3] 组间分离度（各 PC 上 Case 与 Control 的均值差 / 合并 SD）")
    for c in ["PC1", "PC2", "PC3"]:
        a = pc.loc[pc.group == "Case", c]; b = pc.loc[pc.group == "Control", c]
        d = (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        print("    %-4s  Cohen's d = %+.2f   (Case %.1f vs Control %.1f)" % (c, d, a.mean(), b.mean()))

    print("\n  [4] 组织对 PC 的影响（同一 PC 上不同组织的均值）")
    for c in ["PC1", "PC2"]:
        print("    %s: %s" % (c, {k: round(v, 1) for k, v in pc.groupby("tissue")[c].mean().items()}))

    # ---- 离群样本 ----
    print("\n  [5] 离群样本（到组内中心的欧氏距离，>3 倍 MAD 记为离群）")
    out = []
    for c in ["Case", "Control"]:
        sub = pc[pc.group == c]
        cen = sub[["PC1", "PC2", "PC3"]].mean()
        dist = np.sqrt(((sub[["PC1", "PC2", "PC3"]] - cen) ** 2).sum(axis=1))
        med = dist.median(); mad = (dist - med).abs().median() or 1e-9
        for s, dv in zip(sub["sample"], dist):
            z = 0.6745 * (dv - med) / mad
            if z > 3.5:
                out.append((s, c, round(dv, 1), round(z, 2)))
    if out:
        for o in out:
            print("    %-24s %-8s dist=%.1f  robust-z=%.2f" % o)
    else:
        print("    无")

    # ---- 组织匹配后的重分析（用 PCA 前若干成分看差异是否消失）----
    print("\n  [6] 仅比较同一组织内部（Case vs Control）的 PC1 差异")
    for t in pc.tissue.unique():
        sub = pc[pc.tissue == t]
        if sub.group.nunique() < 2 or len(sub) < 6:
            print("    %-14s 样本不足或单组 (n=%d)，跳过" % (t, len(sub)))
            continue
        a = sub.loc[sub.group == "Case", "PC1"]; b = sub.loc[sub.group == "Control", "PC1"]
        d = (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        print("    %-14s n=%2d (Case %2d / Control %2d)  PC1 Cohen's d = %+.2f"
              % (t, len(sub), len(a), len(b), d))
    print()
