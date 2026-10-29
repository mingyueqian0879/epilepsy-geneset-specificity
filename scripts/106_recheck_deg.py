#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 106: 复核手稿中的三个 DEG 率 —— 59.7% / 44.3% / 0.5%
  (a) TLE-HS(64) vs 全部对照(27)                    -> 手稿称 59.7%
  (b) TLE-HS(64) vs ControlCortex(14)，不加组织协变量 -> 手稿称 44.3%
  (c) TLE-HS(64) vs ControlCortex(14)，加组织协变量   -> 手稿称 0.5%
用 limma 与脚本 10/12 相同的设定重跑。
"""
import os
import pandas as pd

BASE = "."
OUT = os.path.join(BASE, "out")
for f in ["expr_matrix.tsv", "design_aligned.csv"]:
    p = os.path.join(OUT, f)
    print("%-22s %s" % (f, "存在 %.1f MB" % (os.path.getsize(p) / 1e6) if os.path.exists(p) else "缺失"))

meta = pd.read_csv(os.path.join(OUT, "design_aligned.csv"))
print("\ndesign_aligned.csv 列:", list(meta.columns))
print("disease:", dict(meta["disease"].value_counts()))
print("tissue :", dict(meta["tissue"].value_counts()))
print("\n按 disease × tissue 交叉:")
print(pd.crosstab(meta["disease"], meta["tissue"]).to_string())
