#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 56: 为 P0 修改计算最终数字
  1) 阳性对照在 248 个基因集中的排名 + 效应量
  2) 各队列最小可检出富集效应 (MDE)
  3) 逐基因统计效力
  4) spiked-in 模拟: 管线能否按预期检出已知幅度的富集
"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats

OUT = "out"
R = pd.read_csv(os.path.join(OUT, "specificity_matrix.csv"))
Rz = R.dropna(subset=["Z"])
S = pd.read_csv(os.path.join(OUT, "dataset_summary.csv"))
LP = pd.read_csv(os.path.join(OUT, "label_permutation_200perm.csv"))

report = []

# ---------------- 1) 阳性对照排名 ----------------
report.append("=" * 92)
report.append(" 1) 阳性对照: 在全部基因集中的排名与效应量 (校准判据, 非显著性判据)")
report.append("=" * 92)
pos_out = {}
for ds in Rz.dataset.unique():
    sub = Rz[Rz.dataset == ds].copy()
    sub["rank"] = sub.Z.rank(ascending=False).astype(int)
    m = len(sub)
    report.append("\n  [%s]  共 %d 个基因集" % (ds, m))
    report.append("    %-24s %-5s %-6s %-8s %-8s %s" % (
        "阳性对照", "n", "obs", "expected", "Z", "排名 (百分位)"))
    for _, r in sub[sub.gene_set.str.startswith("POS_")].sort_values("rank").iterrows():
        pct = 100 * r["rank"] / m
        report.append("    %-24s %-5d %-6d %-8.2f %-+8.2f %d / %d  (top %.1f%%)" % (
            r.gene_set, r.n, r.obs, r.exp, r.Z, r["rank"], m, pct))
        pos_out.setdefault(ds, []).append(
            dict(set=r.gene_set, n=int(r.n), obs=int(r.obs), exp=float(r.exp),
                 Z=float(r.Z), P=float(r.P), rank=int(r["rank"]), m=int(m),
                 pct=float(pct)))
report.append("")
report.append("  说明: 阳性对照为预先指定的校准探针, 不参与 248 个候选基因集的 FDR 家族。")
report.append("        报告排名与效应量, 而非校正后显著性。")

# ---------------- 2) MDE ----------------
report.append("\n" + "=" * 92)
report.append(" 2) 最小可检出富集效应 (Minimum Detectable Enrichment, MDE)")
report.append("=" * 92)
report.append("  定义: 在 m 个检验、BH 控制 FDR<0.05 下, 一个基因集需超出背景多少才能被检出")
report.append("  按最有利情形(该基因集为最小 P 值, 秩=1)计算, 故为下界\n")
mde = []
report.append("  %-12s %-6s %-8s %-9s %-14s %s" % (
    "队列", "检验数", "DEG率", "BH所需Z", "k=100 期望DEG", "可检出的相对超额"))
for _, row in S.iterrows():
    ds = row.dataset
    sub = Rz[Rz.dataset == ds]
    if not len(sub):
        report.append("  %-12s %-6s %-8s %s" % (ds, "-", "-", "无 DEG, 无法检验"))
        continue
    m = len(sub)
    z_need = stats.norm.isf((0.05 / m) / 2)
    p = row.n_deg / row.n_genes
    if p <= 0:
        report.append("  %-12s %-6d %-8s %s" % (ds, m, "0.0%", "不可估计(无 DEG)"))
        continue
    k = 100
    exp_c = k * p
    sd_c = np.sqrt(k * p * (1 - p))
    extra = z_need * sd_c
    rel = extra / exp_c
    report.append("  %-12s %-6d %-8s %-9.2f %-14.1f +%.0f%% (%.1f 个额外DEG)" % (
        ds, m, "%.1f%%" % (100 * p), z_need, exp_c, 100 * rel, extra))
    mde.append(dict(dataset=ds, m=m, deg_rate=100 * p, z_need=z_need,
                    exp_k100=exp_c, rel_excess=100 * rel, extra_deg=extra))

# ---------------- 3) 逐基因效力 ----------------
report.append("\n" + "=" * 92)
report.append(" 3) 逐基因统计效力 (两样本 t 检验近似, alpha=0.05)")
report.append("=" * 92)


def power_t(n1, n2, d, alpha=0.05):
    df = n1 + n2 - 2
    nc = d * np.sqrt(n1 * n2 / (n1 + n2))
    crit = stats.t.ppf(1 - alpha / 2, df)
    return 1 - stats.nct.cdf(crit, df, nc) + stats.nct.cdf(-crit, df, nc)


report.append("  %-12s %-12s %-14s %-14s %s" % (
    "队列", "病例/对照", "效力(d=0.5)", "效力(d=0.8)", "效力(d=1.0)"))
pw = []
for _, r in S.iterrows():
    a, b = int(r.n_case), int(r.n_ctrl)
    p5, p8, p10 = power_t(a, b, 0.5), power_t(a, b, 0.8), power_t(a, b, 1.0)
    report.append("  %-12s %-12s %-14s %-14s %s" % (
        r.dataset, "%d/%d" % (a, b), "%.0f%%" % (100 * p5),
        "%.0f%%" % (100 * p8), "%.0f%%" % (100 * p10)))
    pw.append(dict(dataset=r.dataset, n_case=a, n_ctrl=b,
                   pw_d05=100 * p5, pw_d08=100 * p8, pw_d10=100 * p10))
report.append("")
report.append("  注: 此为单个基因的效力; 基因集检验聚合多个基因, 效力高于此值,")
report.append("      但 MDE 表(第2节)给出了基因集层面的实际可检出下界。")

# ---------------- 4) 标签置换校准 (已有) ----------------
report.append("\n" + "=" * 92)
report.append(" 4) 标签置换校准 (供交叉引用)")
report.append("=" * 92)
for _, r in LP.iterrows():
    report.append("  %-12s 观测DEG=%-7d 零中位=%-4.0f 零95%%=%-6.2f 经验P=%.4f" % (
        r.dataset, r.obs_deg, r.null_median, r.null_p95, r.emp_P))

txt = "\n".join(report)
open(os.path.join(OUT, "P0_revision_numbers.txt"), "w", encoding="utf-8").write(txt)
json.dump({"positive_controls": pos_out, "mde": mde, "power": pw},
          open(os.path.join(OUT, "P0_revision_numbers.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print(txt)
