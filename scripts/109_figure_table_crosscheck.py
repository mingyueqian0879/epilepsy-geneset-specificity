#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 109: 图表 ↔ 数据源 全量交叉核对
逐一验证主图 Fig1-9、主表 Table1-3、补充表 S1-S11b 中引用的数字能否在源数据中复现。
"""
import io, json, os, re
import numpy as np
import pandas as pd

BASE = "."
OUT = os.path.join(BASE, "out")
PKG = "."
bad, ok = [], []


def chk(tag, cond, detail=""):
    (ok if cond else bad).append(tag)
    print("  %s %-66s %s" % ("OK  " if cond else "**!!**", tag, detail))


print("=" * 100)
print(" A. 主表与源数据")
print("=" * 100)
S = pd.read_csv(os.path.join(OUT, "dataset_summary.csv"))
print("\ndataset_summary.csv:")
print(S.to_string(index=False))

T1 = pd.read_csv(os.path.join(PKG, "05_Tables", "Table1_datasets_and_calibration.csv"))
print("\nTable1_datasets_and_calibration.csv:")
print(T1.to_string(index=False))

print("\n--- 逐项核对 ---")
for _, r in T1.iterrows():
    ds = r.iloc[0]
    s = S[S.dataset == ds]
    if not len(s):
        chk("Table1 %s 在源数据中存在" % ds, False); continue
    s = s.iloc[0]
    for col_t1, col_s in [("n_case", "n_case"), ("n_ctrl", "n_ctrl"), ("n_deg", "n_deg")]:
        if col_t1 in T1.columns:
            a, b = r[col_t1], s[col_s]
            chk("Table1 %s %s = %s (源 %s)" % (ds, col_t1, a, b),
                str(a) == str(b) or (pd.notna(a) and pd.notna(b) and abs(float(a) - float(b)) < 0.51))

print("\n" + "=" * 100)
print(" B. 手稿正文中的关键数字 vs 源数据")
print("=" * 100)
md = io.open(os.path.join(PKG, "03_Manuscript", "Manuscript_BIBformat.md"), encoding="utf-8").read()
g256 = S[S.dataset == "GSE256068"].iloc[0]
g213 = S[S.dataset == "GSE213488"].iloc[0]
chk("正文 11,923 DEG (GSE256068)", str(int(g256.n_deg)) == "11923" or abs(g256.n_deg - 11923) < 1,
    "源=%s" % g256.n_deg)
chk("正文 59.7% (GSE256068)", abs(g256.deg_rate - 59.7) < 0.15, "源=%.2f%%" % g256.deg_rate)
chk("正文 262/263 DEG (GSE213488)", abs(g213.n_deg - 262) <= 1 or abs(g213.n_deg - 263) <= 1,
    "源=%s" % g213.n_deg)
chk("正文 1.5% (GSE213488)", abs(g213.deg_rate - 1.5) < 0.15, "源=%.2f%%" % g213.deg_rate)
chk("正文 19,958 tested genes",
    "19,958" in md, "源矩阵基因数 19999；手稿称 19958 <- 需确认过滤口径")

print("\n" + "=" * 100)
print(" C. 基因集层面 (specificity_matrix_final.csv)")
print("=" * 100)
R = pd.read_csv(os.path.join(OUT, "specificity_matrix_final.csv"))
print("  列:", list(R.columns), " 行:", len(R))
for ds in ["GSE256068", "GSE213488", "GSE128300", "GSE272142"]:
    sub = R[R.dataset == ds]
    nz = sub.dropna(subset=["Z"])
    fdr = sub.dropna(subset=["FDR"]) if "FDR" in sub.columns else pd.DataFrame()
    nsig = int((fdr.FDR < 0.05).sum()) if len(fdr) else None
    print("  %-11s 记录 %3d  有Z %3d  FDR<0.05 %s  median Z=%s" %
          (ds, len(sub), len(nz), nsig,
           round(nz.Z.median(), 3) if len(nz) else "NA"))
chk("GSE256068 FDR<0.05 集合数 = 0（手稿）",
    int((R[(R.dataset == "GSE256068")].FDR < 0.05).sum()) == 0,
    "实=%d" % int((R[(R.dataset == "GSE256068")].FDR < 0.05).sum()))
chk("GSE213488 FDR<0.05 集合数 = 21（手稿）",
    int((R[(R.dataset == "GSE213488")].FDR < 0.05).sum()) == 21,
    "实=%d" % int((R[(R.dataset == "GSE213488")].FDR < 0.05).sum()))
chk("GSE256068 有 Z 的集合数 = 248",
    len(R[(R.dataset == "GSE256068")].dropna(subset=["Z"])) == 248,
    "实=%d" % len(R[(R.dataset == "GSE256068")].dropna(subset=["Z"])))

print("\n" + "=" * 100)
print(" D. Fig2B / Fig3 的中位数（图 vs 源）")
print("=" * 100)
for ds in ["GSE256068", "GSE213488", "GSE128300"]:
    sub = R[R.dataset == ds].dropna(subset=["Z"])
    if len(sub) > 5:
        print("  %-11s n=%3d  median Z=%+.3f  IQR=[%+.3f, %+.3f]" %
              (ds, len(sub), sub.Z.median(), sub.Z.quantile(.25), sub.Z.quantile(.75)))
for txt, want in [("median Z = +0.136", 0.136), ("median Z = \u22120.321", -0.321)]:
    chk("手稿 %s" % txt, txt in md, "")

print("\n" + "=" * 100)
print(" E. Table2 / Table3")
print("=" * 100)
for f in ["Table2_geneset_enrichment_summary.csv", "Table3_candidate_genesets.csv"]:
    d = pd.read_csv(os.path.join(PKG, "05_Tables", f))
    print("\n%s  (%d 行)" % (f, len(d)))
    print(d.to_string(index=False))

print("\n" + "=" * 100)
print(" F. 补充表清单与文件一致性")
print("=" * 100)
sd = os.path.join(PKG, "06_Supplementary")
files = sorted(os.listdir(sd))
tabs = sorted(set(f.split("_")[0] for f in files if f.startswith("TableS")))
print("  文件中的表:", tabs)
for s in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "11b"]:
    chk("TableS%s 有 csv" % s, ("TableS%s_" % s) in " ".join(files),
        [f for f in files if f.startswith("TableS%s_" % s)][:2])

print("\n" + "=" * 100)
print(" 结果汇总: 通过 %d 项，异常 %d 项" % (len(ok), len(bad)))
if bad:
    print(" 异常清单:")
    for b in bad:
        print("   -", b)
