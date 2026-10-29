#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step ②-A(重做): 构建 MR 工具变量表 (不做 rsID 映射, 改用 liftover 坐标匹配)
"""
import os, re
import pandas as pd

BASE = "."
OUT  = os.path.join(BASE, "out")
GT   = os.path.join(BASE, "data", "gtex", "GTEx_Analysis_v8_eQTL_independent")

TISSUES = ["Brain_Cortex", "Brain_Hippocampus", "Brain_Frontal_Cortex_BA9",
           "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24"]

lrg_ensg = set()
for line in open(os.path.join(OUT, "LRG_final_ensg.txt"), encoding="utf-8"):
    line = line.strip()
    if line: lrg_ensg.add(line.split(".")[0])
print(f"[LRG] {len(lrg_ensg)} 个 ENSG")

rows = []
for tis in TISSUES:
    fp = os.path.join(GT, f"{tis}.v8.independent_eqtls.txt.gz")
    if not os.path.exists(fp):
        print(f"  [缺失] {tis}"); continue
    d = pd.read_csv(fp, sep="\t")
    d["ensg"] = d["gene_id"].astype(str).str.split(".").str[0]
    hit = d[d["ensg"].isin(lrg_ensg)].copy()
    hit["tissue"] = tis
    rows.append(hit)
    print(f"  {tis:<42} LRG基因={hit['ensg'].nunique():>3}  SNP={len(hit)}")

iv = pd.concat(rows, ignore_index=True)

def parse_vid(v):
    m = re.match(r"chr([\dXYM]+)_(\d+)_([ACGT]+)_([ACGT]+)_b38", str(v))
    return (m.group(1), int(m.group(2)), m.group(3), m.group(4)) if m else (None,None,None,None)

iv[["chrom","pos38","ref","alt"]] = iv["variant_id"].apply(lambda v: pd.Series(parse_vid(v)))
iv = iv[iv["pos38"].notna()].copy()
iv["pos38"] = iv["pos38"].astype(int)

keep = ["ensg","tissue","variant_id","chrom","pos38","ref","alt",
        "tss_distance","maf","pval_nominal","slope","slope_se","pval_beta","rank"]
iv = iv[keep]
iv = iv.sort_values("rank").drop_duplicates(subset=["ensg","tissue","variant_id"])

iv.to_csv(os.path.join(OUT, "MR_instruments.csv"), index=False, encoding="utf-8-sig")
print(f"\n[保存] MR_instruments.csv")
print(f"  行数     = {len(iv)}")
print(f"  基因数   = {iv['ensg'].nunique()}")
print(f"  唯一变异 = {iv['variant_id'].nunique()}")
print(f"  slope 缺失 = {iv['slope'].isna().sum()}   slope_se 缺失 = {iv['slope_se'].isna().sum()}")
print("\n--- 每基因工具变量数分布 ---")
c = iv.groupby("ensg").size()
print(c.describe().to_string())
print("\n--- 每脑区工具变量数 ---")
print(iv.groupby("tissue").size().to_string())
