#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 112: 领域调查、功效表、校准探针表核对"""
import io, json, os, re
import pandas as pd

BASE = "."
OUT = os.path.join(BASE, "out")
PKG = "."
md = io.open(os.path.join(PKG, "03_Manuscript", "Manuscript_BIBformat.md"), encoding="utf-8").read()
body = md[:md.find("## References")]

print("=" * 96)
print(" 1. Fig1 领域调查")
print("=" * 96)
SURV = json.load(io.open(os.path.join(OUT, "field_survey_papers.json"), encoding="utf-8"))
n = len(SURV)
nsp = sum(1 for r in SURV if r.get("has_specificity_test"))
degs = [r["max_deg"] for r in SURV if r.get("max_deg")]
print("  记录数 n            = %d   (手稿称 576)" % n)
print("  有特异性检验         = %d   (手稿称 9, 1.6%%)  -> %.2f%%" % (nsp, 100 * nsp / n))
print("  报告了 DEG 计数的     = %d   (手稿称 70 篇, 12.2%%) -> %.2f%%" % (len(degs), 100 * len(degs) / n))
print("  报告 DEG 最大值 range = %d .. %d   中位数 = %d" % (min(degs), max(degs), int(pd.Series(degs).median())))
print("  手稿称: 12 到 25,101 (中位 326)")
from collections import Counter
c = Counter()
for r in SURV:
    for m in r.get("methods", []):
        c[m] += 1
print("\n  方法频次 top8:")
for k, v in c.most_common(8):
    print("    %-16s %3d  (%.1f%%)" % (k, v, 100 * v / n))

print("\n" + "=" * 96)
print(" 2. TableS7 统计功效")
print("=" * 96)
P7 = pd.read_csv(os.path.join(PKG, "06_Supplementary", "TableS7_statistical_power.csv"))
print(P7.to_string(index=False))
print("  手稿称: d=0.5 时 GSE256068 58%，其余 12-21%")

print("\n" + "=" * 96)
print(" 3. TableS9 校准探针")
print("=" * 96)
P9 = pd.read_csv(os.path.join(PKG, "06_Supplementary", "TableS9_calibration_probes.csv"))
print(P9.to_string(index=False))
print("  手稿称: 小胶质排名 2/248 且 FDR<1e-4 (GSE213488) / 0.234 (GSE256068)")
print("          少突胶质 FDR 0.152 / 0.609；排名 top 12%%")

print("\n" + "=" * 96)
print(" 4. TableS10 spike-in")
print("=" * 96)
P10 = pd.read_csv(os.path.join(PKG, "06_Supplementary", "TableS10_spikein_simulation.csv"))
print(P10.to_string(index=False))

print("\n" + "=" * 96)
print(" 5. TableS11b 方法汇总")
print("=" * 96)
P11 = pd.read_csv(os.path.join(PKG, "06_Supplementary", "TableS11b_method_summary.csv"))
print(P11.to_string(index=False))
print("  手稿称 GSE256068: 0 / 5 / 22 / 148")

print("\n" + "=" * 96)
print(" 6. 手稿中 11,923 / 11,905 / 59.7 的全部出现位置")
print("=" * 96)
for pat in ["11,923", "11,905", "59.7", "59.6"]:
    for m in re.finditer(re.escape(pat), body):
        print("  [%s] ...%s..." % (pat, body[max(0, m.start() - 95):m.start() + 55].replace("\n", " ")))

print("\n" + "=" * 96)
print(" 7. 手稿中 262 / 263 / null 95th 的出现位置")
print("=" * 96)
for pat in ["262", "263", "95th percentile of 20", "percentile of 20"]:
    for m in re.finditer(re.escape(pat), body):
        print("  [%s] ...%s..." % (pat, body[max(0, m.start() - 95):m.start() + 55].replace("\n", " ")))
