#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step ②-C: 两样本 MR 试点 (liftover 版本)
  暴露: GTEx v8 脑 cis-eQTL (hg38) --liftover--> hg19 --> 与 ILAE GWAS 按坐标+等位基因匹配
  结局: ILAE GWAS (beta 列 = Z 值, 需要由 EAF+样本量反推 SE 与 beta)
  方法: lead SNP Wald ratio (主) + 多 SNP IVW (辅)
"""
import gzip, os, re, json, sys
import numpy as np
import pandas as pd
from scipy import stats
from pyliftover import LiftOver

BASE = "."
OUT  = os.path.join(BASE, "out")

GWAS = {
 "Epilepsy_all":   dict(file="GCST90271608.tsv.gz", ncase=29944, nctl=52538, label="癫痫(全部)"),
 "Focal_epilepsy": dict(file="GCST90271610.tsv.gz", ncase=16384, nctl=52538, label="局灶性癫痫"),
 "GGE":            dict(file="GCST90271609.tsv.gz", ncase=7407,  nctl=52538, label="全面性癫痫"),
}
COMP = {"A":"T","T":"A","C":"G","G":"C"}

# ---------- 1. 工具变量 ----------
iv = pd.read_csv(os.path.join(OUT, "MR_instruments.csv"))
def parse_vid(v):
    m = re.match(r"chr([\dXYM]+)_(\d+)_([ACGT]+)_([ACGT]+)_b38", str(v))
    return (m.group(1), int(m.group(2)), m.group(3), m.group(4)) if m else (None,None,None,None)
iv[["chrom","pos38","ref","alt"]] = iv["variant_id"].apply(lambda v: pd.Series(parse_vid(v)))
iv = iv[iv["pos38"].notna()].copy()
iv["rank"] = pd.to_numeric(iv["rank"], errors="coerce")
iv = iv.sort_values(["rank"]).drop_duplicates(subset=["ensg","tissue","variant_id"])
print(f"[工具变量] {len(iv)} 行, {iv['ensg'].nunique()} 基因, {iv['variant_id'].nunique()} 变异")

# ---------- 2. liftover hg38 -> hg19 ----------
chain = os.path.join(BASE, "data", "hg38ToHg19.over.chain.gz")
lo = LiftOver(chain)
uniq_v = iv[["chrom","pos38","ref","alt","variant_id"]].drop_duplicates("variant_id").copy()
print(f"[liftover] {len(uniq_v)} 个变异 hg38 -> hg19 ...")
p19 = []
for c, p in zip(uniq_v["chrom"], uniq_v["pos38"]):
    r = lo.convert_coordinate(f"chr{c}", int(p) - 1)   # pyliftover 用 0-based
    p19.append(int(r[0][1]) + 1 if r else None)        # 转回 1-based
uniq_v["pos19"] = p19
ok = uniq_v["pos19"].notna().sum()
print(f"  成功 {ok} / {len(uniq_v)}  ({100*ok/len(uniq_v):.1f}%)")
vmap = uniq_v.set_index("variant_id")[["pos19"]]
iv = iv.join(vmap, on="variant_id")
iv = iv[iv["pos19"].notna()].copy()
iv["pos19"] = iv["pos19"].astype(int)
iv["key"] = iv["chrom"] + ":" + iv["pos19"].astype(str)
print(f"[对齐后] {len(iv)} 行可用")

# 每个基因的 lead 工具变量(按组织优先级)
TP = ["Brain_Cortex","Brain_Hippocampus","Brain_Frontal_Cortex_BA9",
      "Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24"]
iv["tprio"] = iv["tissue"].map({t:i for i,t in enumerate(TP)}).fillna(9)
lead = iv.sort_values(["tprio","rank"]).drop_duplicates("ensg", keep="first")
print(f"[lead] {len(lead)} 个基因")

key_set = set(iv["key"])

# ---------- 3. 扫描 GWAS ----------
def scan(fname, keys):
    fp = os.path.join(BASE, "data", "gwas", fname)
    out = []
    with gzip.open(fp, "rt", encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {c:i for i,c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr): continue
            k = p[ci["chromosome"]] + ":" + p[ci["base_pair_location"]]
            if k in keys:
                out.append(p)
    d = pd.DataFrame(out, columns=hdr)
    for c in ["beta","effect_allele_frequency","p_value","base_pair_location"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["key"] = d["chromosome"].astype(str) + ":" + d["base_pair_location"].astype(int).astype(str)
    d["effect_allele"] = d["effect_allele"].str.upper()
    d["other_allele"]  = d["other_allele"].str.upper()
    return d

def harmonize(m):
    """把 GWAS 效应统一到 GTEx 的 ALT 等位基因方向"""
    ea, oa, ref, alt = m["effect_allele"], m["other_allele"], m["ref"], m["alt"]
    direct = (ea == alt) & (oa == ref)                 # GWAS 效应等位 = ALT
    flip   = (ea == ref) & (oa == alt)                 # GWAS 效应等位 = REF
    # 链翻转 (互补链)
    ea_c, oa_c = ea.map(COMP), oa.map(COMP)
    sflip_d = (ea_c == alt) & (oa_c == ref)
    sflip_f = (ea_c == ref) & (oa_c == alt)
    keep = direct | flip | sflip_d | sflip_f
    sign = np.where(direct | sflip_d, 1.0, -1.0)
    return keep, sign

allres = {}
for key, cfg in GWAS.items():
    print(f"\n{'='*78}\n[扫描 GWAS] {key}  {cfg['label']}")
    g = scan(cfg["file"], key_set)
    print(f"  坐标匹配到 {len(g)} 个变异")
    if len(g) == 0:
        allres[key] = None; continue

    p_eaf = g["effect_allele_frequency"].clip(0.001, 0.999)
    g["se_out"]   = np.sqrt((1/cfg["ncase"] + 1/cfg["nctl"]) / (2*p_eaf*(1-p_eaf)))
    g["beta_out"] = g["beta"] * g["se_out"]     # beta 列 = Z

    m = iv.merge(g, on="key", how="inner", suffixes=("_iv","_g"))
    print(f"  与工具变量交集: {len(m)} 行")
    if len(m) == 0:
        allres[key] = None; continue

    keep, sign = harmonize(m)
    m = m[keep].copy(); m["sign"] = sign[keep]
    print(f"  等位基因一致: {len(m)} 行  ({m['ensg'].nunique()} 基因)")
    if len(m) == 0:
        allres[key] = None; continue

    m["beta_out_alt"] = m["beta_out"] * m["sign"]
    m["b_exp"] = m["slope"]; m["se_exp"] = m["slope_se"]
    m = m[m["b_exp"].notna() & m["se_exp"].notna() & (m["b_exp"].abs() > 1e-8)].copy()
    m["ratio"]    = m["beta_out_alt"] / m["b_exp"]
    m["ratio_se"] = np.sqrt(m["se_out"]**2/m["b_exp"]**2 +
                            m["beta_out_alt"]**2*m["se_exp"]**2/m["b_exp"]**4)
    m["ratio_z"]  = m["ratio"]/m["ratio_se"]
    m["ratio_p"]  = 2*stats.norm.sf(np.abs(m["ratio_z"]))
    m["F_stat"]   = (m["b_exp"]/m["se_exp"])**2
    m["ci_lo"]    = m["ratio"] - 1.96*m["ratio_se"]
    m["ci_hi"]    = m["ratio"] + 1.96*m["ratio_se"]

    # 基因水平: lead SNP Wald ratio
    gl = m.sort_values(["tprio","rank"]).drop_duplicates("ensg", keep="first").copy()
    gl = gl.sort_values("ratio_p")
    gl["fdr"] = stats.false_discovery_control(gl["ratio_p"].clip(1e-300))

    allres[key] = dict(all=m, gene=gl)
    gl.to_csv(os.path.join(OUT, f"MRpilot_gene_{key}.csv"), index=False, encoding="utf-8-sig")
    m.to_csv(os.path.join(OUT, f"MRpilot_snp_{key}.csv"), index=False, encoding="utf-8-sig")

    print(f"  >>> 完成: 检验 {len(gl)} 个基因 | 名义P<0.05 = {int((gl['ratio_p']<0.05).sum())}"
          f" | FDR<0.05 = {int((gl['fdr']<0.05).sum())}")
    print(f"      F 中位数 = {gl['F_stat'].median():.1f}")

# ---------- 4. 汇总 ----------
print("\n" + "="*78)
print(" MR 试点汇总 (每基因 lead 脑 cis-eQTL, 单 SNP Wald ratio)")
print("="*78)
sym = {}
try:
    mj = json.load(open(os.path.join(OUT, "LRG_map.json"), encoding="utf-8"))
    for s, v in mj.items(): sym[v["ensg"].split(".")[0]] = s
except Exception: pass

for k, r in allres.items():
    if not r: continue
    gl = r["gene"]
    print(f"\n【{GWAS[k]['label']}】 检验 {len(gl)} 基因")
    top = gl[gl["ratio_p"] < 0.05]
    if len(top):
        show = top.head(12).copy()
        show["symbol"] = show["ensg"].map(sym).fillna(show["ensg"])
        cols = ["symbol","tissue","ratio","ci_lo","ci_hi","ratio_p","fdr","F_stat"]
        print(show[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    else:
        print("  无名义显著基因")

json.dump({k: (int(len(v["gene"])) if v else 0) for k,v in allres.items()},
          open(os.path.join(OUT,"MRpilot_summary.json"),"w"), indent=1)
print("\n[保存] MRpilot_gene_*.csv  MRpilot_snp_*.csv  MRpilot_summary.json")
