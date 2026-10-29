#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 28: 检查各数据集格式"""
import gzip, os
import pandas as pd

DL = "data/multi"


def head_gz(fn, n=3, sep=None):
    p = os.path.join(DL, fn)
    print("\n" + "=" * 78)
    print("[%s]  %.2f MB" % (fn, os.path.getsize(p) / 1e6))
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            print("   ", line.rstrip()[:230])
            if i >= n - 1:
                break
    # 行数
    n_row = 0
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
        for _ in fh:
            n_row += 1
    print("    -> 总行数 = %d" % n_row)


head_gz("GSE134697_CountMatrix.tsv.gz")
head_gz("GSE213488_grey_matter_count.csv.gz")
head_gz("GSE128300_RNAstrandedRev_simp_2019.mx.txt.gz")
head_gz("GSE186334_Cx_cyt_CTLvsW0_countTable_edgeR.txt.gz")

print("\n" + "=" * 78)
print("[GSE272142_genes.rawcount.anno.xlsx]")
p = os.path.join(DL, "GSE272142_genes.rawcount.anno.xlsx")
xl = pd.ExcelFile(p)
print("   sheets:", xl.sheet_names)
d = pd.read_excel(p, sheet_name=0, nrows=5)
print("   shape(前5行) =", d.shape)
print("   columns =", list(d.columns)[:12])
print(d.head(3).to_string()[:600])
