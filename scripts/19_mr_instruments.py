#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step ②-A: 准备 MR 输入
 1) 选 LRG 基因在脑组织(Cortex/Hippocampus/FrontalCortex)的 cis-eQTL 工具变量
 2) GTEx variant_id (chr_pos_ref_alt_b38) -> rsID  (Ensembl REST, GRCh38)
 3) 输出 instruments 表
"""
import gzip, os, json, time, re, urllib.request
import pandas as pd

BASE = "."
OUT  = os.path.join(BASE, "out")
GT   = os.path.join(BASE, "data", "gtex", "GTEx_Analysis_v8_eQTL_independent")
HDR  = {"User-Agent": "Mozilla/5.0 (research)", "Content-Type": "application/json"}

TISSUES = ["Brain_Cortex", "Brain_Hippocampus", "Brain_Frontal_Cortex_BA9",
           "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24"]

# LRG ENSG（GTEx gene_id 是 ENSG 带版本）
lrg_ensg = set()
for line in open(os.path.join(OUT, "LRG_final_ensg.txt"), encoding="utf-8"):
    line = line.strip()
    if line:
        lrg_ensg.add(line.split(".")[0])
print(f"[LRG] {len(lrg_ensg)} 个 ENSG")

# ---------- 1. 收集工具变量 ----------
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
    print(f"  {tis:<40} LRG 基因={hit['ensg'].nunique():>3}  SNP={len(hit)}")

iv = pd.concat(rows, ignore_index=True)
iv = iv[["ensg", "tissue", "variant_id", "tss_distance", "maf",
         "pval_nominal", "slope", "slope_se", "pval_beta", "rank"]]
print(f"\n[工具变量合计] {len(iv)} 行,  {iv['ensg'].nunique()} 个基因,  {iv['variant_id'].nunique()} 个唯一 SNP")

# ---------- 2. GTEx variant_id -> rsID (Ensembl REST, GRCh38) ----------
def parse_vid(v):
    m = re.match(r"(chr[\dXYM]+)_(\d+)_([ACGT]+)_([ACGT]+)_b38", v)
    if not m: return None
    return m.group(1).replace("chr", ""), m.group(2), m.group(3), m.group(4)

def ensembl_rsid(chrom, pos):
    u = f"https://rest.ensembl.org/variation/human/{chrom}:{pos}?content-type=application/json"
    for t in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=HDR), timeout=45) as r:
                d = json.load(r)
            return d.get("name") or (d.get("synonyms") or [None])[0]
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(1.5); continue
            return None
        except Exception:
            time.sleep(0.8)
    return None

uniq = sorted(iv["variant_id"].unique())
print(f"\n[Ensembl] 查询 {len(uniq)} 个变异的 rsID ...")
vid2rs = {}
cache_f = os.path.join(OUT, "variant_rsid_cache.json")
cache = json.load(open(cache_f)) if os.path.exists(cache_f) else {}
todo = [v for v in uniq if v not in cache]
print(f"  缓存命中 {len(uniq)-len(todo)},  需查询 {len(todo)}")
for i, v in enumerate(todo, 1):
    p = parse_vid(v)
    rs = ensembl_rsid(p[0], p[1]) if p else None
    cache[v] = rs
    if i % 50 == 0:
        print(f"    {i}/{len(todo)} ...")
        json.dump(cache, open(cache_f, "w"))
    time.sleep(0.08)

json.dump(cache, open(cache_f, "w"))
vid2rs = cache
n_ok = sum(1 for v in uniq if vid2rs.get(v))
print(f"  rsID 获得 {n_ok} / {len(uniq)}")

iv["rsid"] = iv["variant_id"].map(vid2rs)
iv_ok = iv[iv["rsid"].notna()].copy()
iv_ok["chrom"] = iv_ok["variant_id"].map(lambda v: parse_vid(v)[0])
iv_ok["pos38"] = iv_ok["variant_id"].map(lambda v: int(parse_vid(v)[1]))

iv_ok.to_csv(os.path.join(OUT, "MR_instruments.csv"), index=False, encoding="utf-8-sig")
print(f"\n[保存] MR_instruments.csv  ({len(iv_ok)} 行, {iv_ok['ensg'].nunique()} 基因, {iv_ok['rsid'].nunique()} rsID)")
print("\n--- 每个基因的工具变量数 (top 25) ---")
cnt = iv_ok.groupby("ensg").size().sort_values(ascending=False)
print(cnt.head(25).to_string())
