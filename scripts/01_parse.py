#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1: 解析 GSE63808 series matrix
 - 提取样本元数据（phenotype 分组）
 - 提取表达矩阵
 - 保存为 CSV 供 limma 使用
"""
import gzip, os, csv, sys, json
from collections import Counter

BASE = "."
SM   = os.path.join(BASE, "data", "GSE63808_series_matrix.txt.gz")
OUT  = os.path.join(BASE, "out")
os.makedirs(OUT, exist_ok=True)

# ---------- 读取 ----------
meta_lines = []          # 所有 ! 开头的行
table_rows = []          # 表达矩阵
in_table = False
with gzip.open(SM, "rt", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            in_table = False
            continue
        if in_table:
            table_rows.append(line)
        elif line.startswith("!"):
            meta_lines.append(line)

print(f"[读取] 元数据行={len(meta_lines)}  矩阵行={len(table_rows)}")

# ---------- 解析样本特征 ----------
# 同一 tag 可能出现多次 -> 收集为 list
sample_meta = {}   # tag -> list of lists
for line in meta_lines:
    if not line.startswith("!Sample_"):
        continue
    parts = line.split("\t")
    tag = parts[0]
    vals = [p.strip().strip('"') for p in parts[1:]]
    sample_meta.setdefault(tag, []).append(vals)

# 样本数 = Sample_geo_accession 的长度
gsm = sample_meta.get("!Sample_geo_accession", [[]])[0]
n_samples = len(gsm)
print(f"[样本] n = {n_samples}")

# 打印所有 characteristics 通道，寻找分组变量
print("\n[Sample_characteristics 通道]")
char_keys = [k for k in sample_meta if k.startswith("!Sample_characteristics")]
for k in char_keys:
    for i, vals in enumerate(sample_meta[k]):
        cnt = Counter(vals)
        print(f"  {k} [通道{i}]  唯一值={len(cnt)}")
        for v, c in cnt.most_common(20):
            print(f"      {c:>4}  {v}")

# ---------- 确定分组 ----------
def get_char(prefix):
    """从 characteristics 里提取 key: value 形式，返回 list"""
    for k in char_keys:
        for vals in sample_meta[k]:
            if any(v.startswith(prefix) for v in vals):
                return [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    return None

pheno = get_char("phenotype")
if pheno is None:
    print("\n[!] 未找到 phenotype 字段，尝试其它分组字段")
    pheno = get_char("group") or get_char("disease")

if pheno:
    print(f"\n[分组] phenotype 分布: {dict(Counter(pheno))}")
    with open(os.path.join(OUT, "sample_groups.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gsm", "group"])
        for g, p in zip(gsm, pheno):
            w.writerow([g, p])
    print("[保存] sample_groups.csv")

# ---------- 保存表达矩阵 ----------
# 第一行是表头
hdr = table_rows[0].split("\t")
hdr = [h.strip().strip('"') for h in hdr]
print(f"\n[矩阵] 表头列数={len(hdr)}  前3列={hdr[:3]}")

with open(os.path.join(OUT, "expr_matrix.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(hdr)
    for row in table_rows[1:]:
        w.writerow([c.strip().strip('"') for c in row.split("\t")])
print(f"[保存] expr_matrix.csv  ({len(table_rows)-1} 行探针)")

# ---------- 保存完整元数据 ----------
with open(os.path.join(OUT, "sample_meta.json"), "w", encoding="utf-8") as f:
    json.dump({k: v for k, v in sample_meta.items()}, f, ensure_ascii=False, indent=1)
print("[保存] sample_meta.json")
