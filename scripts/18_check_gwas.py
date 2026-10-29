#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""检查 ILAE GWAS 三份文件的数据完整性"""
import pandas as pd, numpy as np

for g in ["GCST90271608", "GCST90271609", "GCST90271610"]:
    f = rf"data/gwas/{g}.tsv.gz"
    d = pd.read_csv(f, sep="\t", nrows=400000)
    se_na = d["standard_error"].isna().mean()
    print(f"\n{g}:  读 {len(d)} 行")
    print(f"   SE 缺失率   = {se_na:.3f}")
    print(f"   beta 范围   = [{d['beta'].min():.3f}, {d['beta'].max():.3f}]")
    print(f"   EAF 范围    = [{d['effect_allele_frequency'].min():.3f}, {d['effect_allele_frequency'].max():.3f}]")
    print(f"   p 值范围    = [{d['p_value'].min():.2e}, {d['p_value'].max():.3f}]")
    print(f"   等位基因    = {sorted(set(d['effect_allele'].dropna().unique()))[:8]}")
    rs = d["rs_id"].dropna()
    print(f"   rs_id 非空率= {len(rs)/len(d):.3f}  示例={rs.iloc[0] if len(rs) else 'NA'}")
    # 用 beta/p 反推 SE 的可行性
    ok = d["standard_error"].notna() & d["beta"].notna() & d["p_value"].notna()
    if ok.sum() > 100:
        from scipy import stats
        z = np.abs(d.loc[ok, "beta"]) / d.loc[ok, "standard_error"]
        p_implied = 2 * stats.norm.sf(z)
        print(f"   [验证] SE 与 p 一致性: 相关系数 = {np.corrcoef(-np.log10(p_implied.clip(1e-300)), -np.log10(d.loc[ok,'p_value']))[0,1]:.4f}")
    # 若 SE 缺失, 用 beta/p 反推
    miss = d["standard_error"].isna() & d["beta"].notna() & d["p_value"].notna() & (d["p_value"] > 0)
    if miss.sum() > 0:
        from scipy import stats
        z2 = np.abs(stats.norm.isf(d.loc[miss, "p_value"] / 2))
        se_imp = (np.abs(d.loc[miss, "beta"]) / z2)
        print(f"   [可反推] SE 缺失但 beta/p 齐全的行 = {miss.sum()}  反推 SE 中位数 = {se_imp.median():.4f}")
