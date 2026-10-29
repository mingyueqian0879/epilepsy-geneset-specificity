#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 31: 从 Gencode GTF 建 ENSG -> symbol 全基因组映射"""
import gzip, os, re, json
import pandas as pd

BASE = "."
GTF = os.path.join(BASE, "data", "gencode.v44.annotation.gtf.gz")
OUT = os.path.join(BASE, "out", "multi")
os.makedirs(OUT, exist_ok=True)

pat_id = re.compile(r'gene_id "([^"]+)"')
pat_nm = re.compile(r'gene_name "([^"]+)"')
pat_tp = re.compile(r'gene_type "([^"]+)"')

rows = []
with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        attr = p[8]
        gid = pat_id.search(attr)
        gnm = pat_nm.search(attr)
        gtp = pat_tp.search(attr)
        if gid and gnm:
            rows.append((gid.group(1).split(".")[0], gnm.group(1),
                         gtp.group(1) if gtp else ""))

df = pd.DataFrame(rows, columns=["ensg", "symbol", "type"])
df = df.drop_duplicates("ensg")
print("[Gencode] 基因数 =", len(df))
print(df["type"].value_counts().head(8).to_string())

df.to_csv(os.path.join(OUT, "ensg2symbol_gencode.csv"), index=False)
json.dump(dict(zip(df["ensg"], df["symbol"])),
          open(os.path.join(OUT, "ensg2symbol_gencode.json"), "w"), ensure_ascii=False)
print("[saved] ensg2symbol_gencode.csv / .json")

# 检查 4 个数据集的覆盖率
print("\n[覆盖率检查]")
for ds in ["GSE256068", "GSE213488", "GSE128300", "GSE272142"]:
    f = os.path.join(OUT, ds + "_expr.pkl")
    if not os.path.exists(f):
        continue
    e = pd.read_pickle(f)
    ids = set(str(i).split(".")[0] for i in e.index)
    hit = len(ids & set(df["ensg"]))
    print("  %-12s 基因=%-7d 映射命中=%-7d (%.1f%%)" % (ds, len(ids), hit, 100 * hit / len(ids)))
