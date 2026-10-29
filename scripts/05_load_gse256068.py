#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 2: 解析 GSE256068
 - 从 series matrix 取元数据（Sample_description = 矩阵列名）
 - 解析 processed 矩阵（空格分隔、带引号、log2 TPM）
 - 输出 design 表 + 表达矩阵（基因×样本）
"""
import gzip, os, io, csv, json, re
import numpy as np
import pandas as pd

BASE = "."
SM   = os.path.join(BASE, "data", "GSE256068_series_matrix.txt.gz")
PROC = os.path.join(BASE, "data", "GSE256068_processed_data.csv.gz")
OUT  = os.path.join(BASE, "out")
os.makedirs(OUT, exist_ok=True)

# ---------------- 1. 元数据 ----------------
tags = {}
with gzip.open(SM, "rt", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        if not line.startswith("!Sample_"):
            continue
        parts = line.rstrip("\n").split("\t")
        tag = parts[0]
        vals = [p.strip().strip('"') for p in parts[1:]]
        tags.setdefault(tag, []).append(vals)

n = len(tags["!Sample_geo_accession"][0])
print(f"[样本数] {n}")

def col(prefix):
    """合并同一 prefix 的所有通道（取第一个含该 prefix 的通道）"""
    for k, chans in tags.items():
        if k.startswith(prefix):
            for v in chans:
                if v and any(str(x).startswith(prefix.split("_ch")[0]) or True for x in v[:1]):
                    if len(v) == n:
                        return v
    return None

def get_char(key):
    for k, chans in tags.items():
        if "characteristics" not in k:
            continue
        for v in chans:
            if len(v) == n and v[0].lower().startswith(key.lower()):
                return [x.split(":", 1)[1].strip() if ":" in x else x for x in v]
    return None

meta = pd.DataFrame({
    "gsm":        tags["!Sample_geo_accession"][0],
    "title":      tags["!Sample_title"][0],
    "matrix_col": tags["!Sample_description"][0],   # 与 processed 矩阵列名一致
    "tissue":     get_char("tissue"),
    "sex":        get_char("Sex"),
    "age":        get_char("age"),
    "disease":    get_char("disease state"),
})
print("\n[分组统计]")
print(meta.groupby(["disease", "tissue"]).size().to_string())
meta.to_csv(os.path.join(OUT, "design_GSE256068.csv"), index=False, encoding="utf-8-sig")
print(f"\n[保存] design_GSE256068.csv")

# ---------------- 2. 表达矩阵 ----------------
print("\n[解析表达矩阵] ...")
rows = []
with gzip.open(PROC, "rt", encoding="utf-8", errors="replace") as fh:
    header = None
    for i, line in enumerate(fh):
        # 格式: "GENE" v1 v2 ... v162  (引号仅用于标识符)
        parts = line.replace('"', " ").split()
        if not parts:
            continue
        if i == 0:
            header = parts                      # 162 个样本名（无基因列表头）
            print(f"  表头列数 = {len(header)}   前3 = {header[:3]}")
        else:
            rows.append(parts)
        if i % 20000 == 0 and i > 0:
            print(f"  ... {i} 行")

genes = [r[0] for r in rows]
mat = np.array([[float(x) for x in r[1:]] for r in rows], dtype=np.float32)
print(f"  矩阵形状 = {mat.shape}  (基因 × 样本)")

expr = pd.DataFrame(mat, index=genes, columns=header)
# 对齐：只保留 design 中存在的列
common = [c for c in meta["matrix_col"] if c in expr.columns]
print(f"  与 design 匹配的列 = {len(common)} / {len(meta)}")
expr = expr[common]
meta = meta.set_index("matrix_col").loc[common].reset_index()

expr.to_pickle(os.path.join(OUT, "expr_log2tpm.pkl"))
meta.to_csv(os.path.join(OUT, "design_aligned.csv"), index=False, encoding="utf-8-sig")
print(f"\n[保存] expr_log2tpm.pkl  {expr.shape}")
print(f"[保存] design_aligned.csv")

print("\n[表达值范围] min=%.2f  median=%.2f  max=%.2f" %
      (np.nanmin(mat), np.nanmedian(mat), np.nanmax(mat)))
print("[矩阵列名示例]", list(expr.columns[:3]))
print("[design 对应]")
print(meta[["gsm", "title", "disease", "tissue", "age"]].head(6).to_string(index=False))
