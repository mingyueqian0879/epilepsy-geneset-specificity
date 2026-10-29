#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step ②-D: MR 背景对照 (决定性)
  对【全部】有脑 eQTL 的基因跑同样的 MR, 比较 LRG vs 非LRG 的命中率
  -> 直接回答: 乳酸化基因的 MR 信号是否有特异性?
"""
import gzip, os, re, json
import numpy as np
import pandas as pd
from scipy import stats
from pyliftover import LiftOver

BASE = "."
OUT  = os.path.join(BASE, "out")
GT   = os.path.join(BASE, "data", "gtex", "GTEx_Analysis_v8_eQTL_independent")
COMP = {"A":"T","T":"A","C":"G","G":"C"}

TISSUES = ["Brain_Cortex","Brain_Hippocampus","Brain_Frontal_Cortex_BA9",
           "Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24"]
GWAS = {
 "Epilepsy_all":   dict(file="GCST90271608.tsv.gz", ncase=29944, nctl=52538, label="癫痫(全部)"),
 "Focal_epilepsy": dict(file="GCST90271610.tsv.gz", ncase=16384, nctl=52538, label="局灶性癫痫"),
 "GGE":            dict(file="GCST90271609.tsv.gz", ncase=7407,  nctl=52538, label="全面性癫痫"),
}

lrg_ensg = set(l.split(".")[0].strip() for l in open(os.path.join(OUT,"LRG_final_ensg.txt"),encoding="utf-8") if l.strip())
print(f"[LRG] {len(lrg_ensg)} 个 ENSG")

def parse_vid(v):
    m = re.match(r"chr([\dXYM]+)_(\d+)_([ACGT]+)_([ACGT]+)_b38", str(v))
    return (m.group(1), int(m.group(2)), m.group(3), m.group(4)) if m else (None,None,None,None)

# ---------- 1. 全部脑 eQTL (不只 LRG) ----------
rows = []
for tis in TISSUES:
    fp = os.path.join(GT, f"{tis}.v8.independent_eqtls.txt.gz")
    d = pd.read_csv(fp, sep="\t")
    d["ensg"] = d["gene_id"].astype(str).str.split(".").str[0]
    d["tissue"] = tis
    rows.append(d[["ensg","tissue","variant_id","slope","slope_se","rank","maf","pval_nominal"]])
alliv = pd.concat(rows, ignore_index=True)
alliv = alliv.dropna(subset=["slope","slope_se"])
alliv = alliv[alliv["slope"].abs() > 1e-8]
print(f"[全部工具变量] {len(alliv)} 行, {alliv['ensg'].nunique()} 基因, {alliv['variant_id'].nunique()} 变异")

alliv[["chrom","pos38","ref","alt"]] = alliv["variant_id"].apply(lambda v: pd.Series(parse_vid(v)))
alliv = alliv[alliv["pos38"].notna()].copy()

# ---------- 2. liftover ----------
lo = LiftOver(os.path.join(BASE,"data","hg38ToHg19.over.chain.gz"))
uv = alliv[["chrom","pos38","variant_id"]].drop_duplicates("variant_id")
print(f"[liftover] {len(uv)} 个唯一变异 ...")
uv["pos19"] = [ (lambda r: int(r[0][1])+1 if r else None)(lo.convert_coordinate(f"chr{c}", int(p)-1))
                for c,p in zip(uv["chrom"], uv["pos38"]) ]
print(f"  成功 {uv['pos19'].notna().sum()} / {len(uv)}")
alliv = alliv.merge(uv[["variant_id","pos19"]], on="variant_id", how="left")
alliv = alliv[alliv["pos19"].notna()].copy()
alliv["pos19"] = alliv["pos19"].astype(int)
alliv["key"] = alliv["chrom"] + ":" + alliv["pos19"].astype(str)
keyset = set(alliv["key"])
print(f"[可用] {len(alliv)} 行, {alliv['ensg'].nunique()} 基因")

def scan(fname, keys):
    out = []
    with gzip.open(os.path.join(BASE,"data","gwas",fname), "rt", encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {c:i for i,c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr): continue
            if p[ci["chromosome"]] + ":" + p[ci["base_pair_location"]] in keys:
                out.append(p)
    d = pd.DataFrame(out, columns=hdr)
    for c in ["beta","effect_allele_frequency","base_pair_location"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["key"] = d["chromosome"].astype(str)+":"+d["base_pair_location"].astype(int).astype(str)
    d["effect_allele"] = d["effect_allele"].str.upper()
    d["other_allele"]  = d["other_allele"].str.upper()
    return d

summary = []
for gk, cfg in GWAS.items():
    print(f"\n{'='*78}\n[{cfg['label']}] 扫描 GWAS ...")
    g = scan(cfg["file"], keyset)
    print(f"  坐标匹配 {len(g)} 变异")
    p_eaf = g["effect_allele_frequency"].clip(0.001,0.999)
    g["se_out"] = np.sqrt((1/cfg["ncase"]+1/cfg["nctl"])/(2*p_eaf*(1-p_eaf)))
    g["beta_out"] = g["beta"]*g["se_out"]

    m = alliv.merge(g, on="key", how="inner", suffixes=("_iv","_g"))
    ea,oa,ref,alt = m["effect_allele"],m["other_allele"],m["ref"],m["alt"]
    direct=(ea==alt)&(oa==ref); flip=(ea==ref)&(oa==alt)
    eac,oac = ea.map(COMP), oa.map(COMP)
    sf_d=(eac==alt)&(oac==ref); sf_f=(eac==ref)&(oac==alt)
    keep = direct|flip|sf_d|sf_f
    m = m[keep].copy()
    m["sign"] = np.where((direct|sf_d)[keep], 1.0, -1.0)
    m["beta_out_alt"] = m["beta_out"]*m["sign"]
    m["b_exp"]=m["slope"]; m["se_exp"]=m["slope_se"]
    m["ratio"]=m["beta_out_alt"]/m["b_exp"]
    m["ratio_se"]=np.sqrt(m["se_out"]**2/m["b_exp"]**2 + m["beta_out_alt"]**2*m["se_exp"]**2/m["b_exp"]**4)
    m["ratio_z"]=m["ratio"]/m["ratio_se"]
    m["ratio_p"]=2*stats.norm.sf(np.abs(m["ratio_z"]))
    m["F_stat"]=(m["b_exp"]/m["se_exp"])**2
    m["is_lrg"]=m["ensg"].isin(lrg_ensg)

    gl = m.sort_values("rank").drop_duplicates("ensg", keep="first").copy()
    gl = gl[gl["F_stat"] > 10]                      # 仅保留强工具变量
    print(f"  基因水平可检验: {len(gl)}  (LRG {int(gl['is_lrg'].sum())} / 非LRG {int((~gl['is_lrg']).sum())})")

    a = int(((gl["ratio_p"]<0.05) & gl["is_lrg"]).sum());  b = int(gl["is_lrg"].sum())-a
    c_ = int(((gl["ratio_p"]<0.05) & ~gl["is_lrg"]).sum()); d_ = int((~gl["is_lrg"]).sum())-c_
    ft = stats.fisher_exact([[a,b],[c_,d_]], alternative="greater")
    r_lrg = a/max(gl["is_lrg"].sum(),1); r_bg = c_/max((~gl["is_lrg"]).sum(),1)
    print(f"  LRG  命中 {a}/{int(gl['is_lrg'].sum())} = {100*r_lrg:.1f}%")
    print(f"  非LRG 命中 {c_}/{int((~gl['is_lrg']).sum())} = {100*r_bg:.1f}%")
    print(f"  >>> Fisher OR = {ft.statistic:.3f}   P = {ft.pvalue:.3f}   "
          f"{'✅ LRG 富集' if ft.pvalue<0.05 else '❌ 无富集'}")
    # 效应量分布比较
    if gl["is_lrg"].sum()>3:
        u = stats.mannwhitneyu(np.abs(gl.loc[gl['is_lrg'],'ratio_z']),
                               np.abs(gl.loc[~gl['is_lrg'],'ratio_z']))
        print(f"  |Z| 分布 LRG 中位={np.abs(gl.loc[gl['is_lrg'],'ratio_z']).median():.3f} "
              f"vs 非LRG={np.abs(gl.loc[~gl['is_lrg'],'ratio_z']).median():.3f}  (MWU P={u.pvalue:.3f})")
    summary.append(dict(gwas=cfg["label"], n_tested=len(gl), lrg_tested=int(gl['is_lrg'].sum()),
                        lrg_hit=a, bg_hit=c_, rate_lrg=100*r_lrg, rate_bg=100*r_bg,
                        fisher_or=float(ft.statistic), fisher_p=float(ft.pvalue)))
    gl.to_csv(os.path.join(OUT, f"MRbg_gene_{gk}.csv"), index=False, encoding="utf-8-sig")

pd.DataFrame(summary).to_csv(os.path.join(OUT,"MR_background_test.csv"), index=False, encoding="utf-8-sig")
print("\n" + "="*78)
print(" 背景对照汇总: LRG vs 全部其它基因 的 MR 命中率")
print("="*78)
print(pd.DataFrame(summary).to_string(index=False, float_format=lambda x: f"{x:.4g}"))
