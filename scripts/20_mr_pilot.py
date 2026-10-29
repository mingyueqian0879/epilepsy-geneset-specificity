#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step ②-B: 两样本 MR 试点  —— LRG 脑 eQTL 工具变量 vs ILAE 癫痫 GWAS
  暴露: GTEx v8 脑组织 cis-eQTL (slope, slope_se)   [GRCh38, alt 等位基因]
  结局: ILAE GWAS (beta 列 = Z 值, EAF)             [GRCh37, 用 rsID 匹配]
  方法: 单 SNP Wald ratio (主力) + 多 SNP IVW (稳健性)
  关键: 结局 SE 由 EAF 与样本量反推 (SE = sqrt((1/Ncase+1/Nctl)/(2*p*(1-p))))
"""
import gzip, os, re, json
import numpy as np
import pandas as pd
from scipy import stats

BASE = "."
OUT  = os.path.join(BASE, "out")
os.makedirs(OUT, exist_ok=True)

# ---------- 结局样本量 ----------
GWAS = {
 "Epilepsy_all":   dict(file="GCST90271608.tsv.gz", ncase=29944, nctl=52538, label="癫痫(全部)"),
 "Focal_epilepsy": dict(file="GCST90271610.tsv.gz", ncase=16384, nctl=52538, label="局灶性癫痫"),
 "GGE":            dict(file="GCST90271609.tsv.gz", ncase=7407,  nctl=52538, label="全面性癫痫"),
}

# ---------- 1. 读工具变量 ----------
iv = pd.read_csv(os.path.join(OUT, "MR_instruments.csv"))
print(f"[工具变量] {len(iv)} 行, {iv['ensg'].nunique()} 基因, {iv['rsid'].nunique()} rsID")

def parse_vid(v):
    m = re.match(r"chr([\dXYM]+)_(\d+)_([ACGT]+)_([ACGT]+)_b38", str(v))
    return (m.group(1), int(m.group(2)), m.group(3), m.group(4)) if m else (None, None, None, None)

iv[["chrom","pos","ref","alt"]] = iv["variant_id"].apply(lambda v: pd.Series(parse_vid(v)))
iv = iv[iv["rsid"].notna() & iv["alt"].notna()].copy()
iv["rsid"] = iv["rsid"].astype(str)

# 每个基因-组织保留最强(rank 最小)的独立 eQTL
iv["rank"] = pd.to_numeric(iv["rank"], errors="coerce")
iv = iv.sort_values("rank").drop_duplicates(subset=["ensg","tissue","rsid"])

TISSUE_PRIORITY = ["Brain_Cortex", "Brain_Hippocampus", "Brain_Frontal_Cortex_BA9",
                   "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24"]
iv["tprio"] = iv["tissue"].map({t:i for i,t in enumerate(TISSUE_PRIORITY)}).fillna(9)

# 基因水平最佳工具变量 (优先皮层/海马)
lead = (iv.sort_values(["tprio","rank"])
          .drop_duplicates(subset=["ensg"], keep="first")
          .copy())
print(f"[lead 工具变量] {len(lead)} 个基因 (每基因 1 个 lead cis-eQTL)")

need_rs = set(lead["rsid"]) | set(iv["rsid"])
print(f"[GWAS 查找] 需要匹配 {len(need_rs)} 个 rsID")

# ---------- 2. 流式扫描 GWAS ----------
def scan_gwas(fname, rs_set):
    rows = []
    fp = os.path.join(BASE, "data", "gwas", fname)
    with gzip.open(fp, "rt", encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c:i for i,c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr): continue
            rs = p[idx["rs_id"]]
            if rs in rs_set:
                rows.append(p)
    d = pd.DataFrame(rows, columns=hdr)
    for c in ["beta","standard_error","effect_allele_frequency","p_value","base_pair_location"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d

results = {}
for key, cfg in GWAS.items():
    print(f"\n[扫描 GWAS] {key} ({cfg['label']}) ...")
    g = scan_gwas(cfg["file"], need_rs)
    print(f"  匹配到 {len(g)} 个变异")
    g["effect_allele"] = g["effect_allele"].str.upper()
    g["other_allele"]  = g["other_allele"].str.upper()

    # 结局 SE 反推 (Z -> beta)
    p_eaf = g["effect_allele_frequency"].clip(0.001, 0.999)
    se_out = np.sqrt((1/cfg["ncase"] + 1/cfg["nctl"]) / (2 * p_eaf * (1 - p_eaf)))
    g["se_out"]   = se_out
    g["beta_out"] = g["beta"] * se_out        # beta 列是 Z

    # ---------- 3. Harmonize + Wald ratio ----------
    m = lead.merge(g, on="rsid", how="inner", suffixes=("_iv","_out"))
    print(f"  与 lead 工具变量交集: {len(m)} 个基因")
    if len(m) == 0:
        results[key] = None; continue

    # 等位基因对齐: GTEx slope 是 alt 等位基因效应
    # 若 GWAS effect_allele == iv.alt -> 同向; == iv.ref -> 反向; 否则丢弃
    same = m["effect_allele"] == m["alt"]
    flip = m["effect_allele"] == m["ref"]
    m = m[same | flip].copy()
    m["sign"] = np.where(m["effect_allele"] == m["alt"], 1.0, -1.0)
    # 统一到 alt 等位基因方向
    m["beta_out_alt"] = m["beta_out"] * m["sign"]
    print(f"  等位基因对齐后: {len(m)} 个基因  (丢弃 {int((~(same|flip)).sum())})")

    if len(m) == 0:
        results[key] = None; continue

    # Wald ratio
    m["b_exp"] = m["slope"]; m["se_exp"] = m["slope_se"]
    m = m[m["b_exp"].notna() & (m["b_exp"].abs() > 1e-8) & m["se_exp"].notna()]
    m["ratio"]     = m["beta_out_alt"] / m["b_exp"]
    m["ratio_se"]  = np.sqrt(m["se_out"]**2 / m["b_exp"]**2
                             + m["beta_out_alt"]**2 * m["se_exp"]**2 / m["b_exp"]**4)
    m["ratio_z"]   = m["ratio"] / m["ratio_se"]
    m["ratio_p"]   = 2 * stats.norm.sf(np.abs(m["ratio_z"]))
    m["F_stat"]    = (m["b_exp"] / m["se_exp"])**2          # 工具变量强度
    m["ci_lo"]     = m["ratio"] - 1.96 * m["ratio_se"]
    m["ci_hi"]     = m["ratio"] + 1.96 * m["ratio_se"]

    m = m.sort_values("ratio_p")
    m["fdr"] = np.nan
    ok = m["ratio_p"].notna()
    if ok.sum() > 0:
        m.loc[ok, "fdr"] = stats.false_discovery_control(m.loc[ok, "ratio_p"])

    results[key] = m
    m.to_csv(os.path.join(OUT, f"MR_{key}.csv"), index=False, encoding="utf-8-sig")

    n_sig = int((m["ratio_p"] < 0.05).sum())
    n_fdr = int((m["fdr"] < 0.05).sum())
    print(f"  >>> MR 完成: 检验 {len(m)} 个基因 | 名义 P<0.05 = {n_sig} | FDR<0.05 = {n_fdr}")
    print(f"      工具变量强度 F 中位数 = {m['F_stat'].median():.1f}  (F>10 为强工具)")

# ---------- 4. 汇总 ----------
print("\n" + "="*78)
print(" MR 试点结果汇总 (lead cis-eQTL, 单 SNP Wald ratio)")
print("="*78)
for key, m in results.items():
    if m is None: continue
    sig = m[m["ratio_p"] < 0.05].head(10)
    print(f"\n【{GWAS[key]['label']}】 名义显著 {int((m['ratio_p']<0.05).sum())} / {len(m)}")
    if len(sig):
        cols = ["ensg","tissue","rsid","b_exp","ratio","ci_lo","ci_hi","ratio_p","fdr","F_stat"]
        print(sig[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
