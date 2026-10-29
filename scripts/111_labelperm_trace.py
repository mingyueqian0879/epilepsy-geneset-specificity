#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 111: label-permutation 数字溯源（手稿 vs Table1 vs 各源文件）"""
import io, os
import pandas as pd

BASE = "."
OUT = os.path.join(BASE, "out")
PKG = "."

print("=" * 100)
print(" 所有 label permutation 相关文件")
print("=" * 100)
for f in sorted(os.listdir(OUT)):
    if "label" in f.lower() or "perm" in f.lower():
        p = os.path.join(OUT, f)
        try:
            d = pd.read_csv(p)
            print("\n--- out/%s  (%d 行) ---" % (f, len(d)))
            print(d.to_string(index=False))
        except Exception as e:
            print("\n--- out/%s  (读取失败 %s) ---" % (f, str(e)[:40]))

print("\n" + "=" * 100)
print(" 投稿包中的 S5 / S6")
print("=" * 100)
for f in ["TableS5_label_permutation.csv", "TableS6_label_permutation_robustness.csv"]:
    p = os.path.join(PKG, "06_Supplementary", f)
    d = pd.read_csv(p)
    print("\n--- %s ---" % f)
    print(d.to_string(index=False))

print("\n" + "=" * 100)
print(" Table 1 中的 label-permutation 列")
print("=" * 100)
T1 = pd.read_csv(os.path.join(PKG, "05_Tables", "Table1_datasets_and_calibration.csv"))
cols = [c for c in T1.columns if "Null" in c or "Empirical" in c]
print(T1[["Dataset"] + cols].to_string(index=False))

print("\n" + "=" * 100)
print(" Fig5 绘图脚本引用的是哪个文件")
print("=" * 100)
import re
for s in ["39_fig5_redo.py", "38_figures2.py", "66_figs_update.py"]:
    p = os.path.join(BASE, "scripts", s)
    if not os.path.exists(p):
        continue
    t = io.open(p, encoding="utf-8").read()
    hits = re.findall(r'(label_permutation[\w.]*|TableS\d[\w.]*|null_deg_p95|null_p95|p_deg|empP|emp_P)', t)
    print("  %-22s 引用: %s" % (s, sorted(set(hits)) or "无"))
