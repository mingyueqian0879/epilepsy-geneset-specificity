#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step ①: LRG 基因在 GTEx 脑组织中的 cis-eQTL 工具变量覆盖率
"""
import gzip, os, json, glob
import pandas as pd

BASE = "."
GtexDir = os.path.join(BASE, "data", "gtex", "GTEx_Analysis_v8_eQTL_independent")
if not os.path.isdir(GtexDir):
    for c in glob.glob(os.path.join(BASE, "data", "gtex", "**", "*.independent_eqtls.txt.gz"), recursive=True):
        GtexDir = os.path.dirname(c); break
OUT = os.path.join(BASE, "out")
print("[GTEx 目录]", GtexDir)

## ---------- 1. 看格式 ----------
f = os.path.join(GtexDir, "Brain_Hippocampus.v8.independent_eqtls.txt.gz")
df = pd.read_csv(f, sep="\t", nrows=5)
print("\n[文件格式]", os.path.basename(f))
print("  列:", list(df.columns))
print(df.head(3).to_string())

## ---------- 2. 全脑区扫描 ----------
BRAIN = ["Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus",
         "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24",
         "Brain_Caudate_basal_ganglia", "Brain_Cerebellum",
         "Brain_Cerebellar_Hemisphere", "Brain_Hypothalamus",
         "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
         "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra"]

# LRG
lrg = [l.strip() for l in open(os.path.join(OUT, "LRG_final_symbols.txt"), encoding="utf-8") if l.strip()]
print(f"\n[LRG] {len(lrg)} 个基因")

per_tissue = {}
gene_hits = {}          # symbol -> {tissue: n_snp}
all_rows = []
for tis in BRAIN:
    fp = os.path.join(GtexDir, f"{tis}.v8.independent_eqtls.txt.gz")
    if not os.path.exists(fp):
        print(f"  [缺失] {tis}"); continue
    d = pd.read_csv(fp, sep="\t")
    # gene 列名探测
    gcol = None
    for cand in ["gene_id", "gene", "gene_symbol", "phenotype_id"]:
        if cand in d.columns:
            gcol = cand; break
    if gcol is None:
        print(f"  [无法识别基因列] {tis}: {list(d.columns)}"); continue
    d["symbol"] = d[gcol].astype(str).str.split(".").str[0]
    # 若 gene_id 是 ENSG，需要用 symbol 匹配 -> 用 LRG ENSG 反查
    is_ensg = d["symbol"].str.startswith("ENSG").any()
    hit = d[d["symbol"].isin(lrg)] if not is_ensg else None
    if is_ensg:
        lrg_ensg = set(x.strip().split(".")[0] for x in
                       open(os.path.join(OUT, "LRG_final_ensg.txt"), encoding="utf-8") if x.strip())
        hit = d[d["symbol"].isin(lrg_ensg)]
    per_tissue[tis] = {"n_total_genes": d["symbol"].nunique(), "n_hit": hit["symbol"].nunique(),
                       "n_rows_hit": len(hit), "is_ensg": bool(is_ensg)}
    for s, grp in hit.groupby("symbol"):
        gene_hits.setdefault(s, {})[tis] = len(grp)
    all_rows.append(hit.assign(tissue=tis))
    print(f"  {tis:<40} 基因={d['symbol'].nunique():>6}  LRG命中={hit['symbol'].nunique():>4}  SNPs={len(hit)}")

print(f"\n{'='*74}")
print(" 问题① 结论: LRG 在脑组织的 cis-eQTL 工具变量覆盖")
print(f"{'='*74}")
covered = set(gene_hits.keys())
print(f"  LRG 总数           : {len(lrg)}")
print(f"  **有脑 eQTL 的 LRG : {len(covered)}   ({100*len(covered)/len(lrg):.1f}%)**")
multi = {g: t for g, t in gene_hits.items() if len(t) >= 2}
print(f"  在 >=2 个脑区有 eQTL: {len(multi)}")
core = {g: t for g, t in gene_hits.items() if {"Brain_Cortex"} <= set(t) or {"Brain_Hippocampus"} <= set(t)}
print(f"  在 Cortex 或 Hippocampus 有 eQTL: {len(core)}")

# 汇总表
summ = []
for g, t in gene_hits.items():
    summ.append({"symbol": g, "n_tissues": len(t),
                 "total_snps": sum(t.values()),
                 "tissues": ",".join(sorted(k.replace("Brain_", "") for k in t))})
S = pd.DataFrame(summ).sort_values(["total_snps"], ascending=False)
S.to_csv(os.path.join(OUT, "LRG_brain_eqtl_coverage.csv"), index=False, encoding="utf-8-sig")
print(f"\n[保存] LRG_brain_eqtl_coverage.csv")
print("\n--- 工具变量最多的 25 个 LRG ---")
print(S.head(25).to_string(index=False))

if all_rows:
    pd.concat(all_rows).to_csv(os.path.join(OUT, "LRG_brain_eqtl_all.csv"), index=False)
    print(f"\n[保存] LRG_brain_eqtl_all.csv")
json.dump(per_tissue, open(os.path.join(OUT, "gtex_coverage_summary.json"), "w"), indent=1, ensure_ascii=False)
