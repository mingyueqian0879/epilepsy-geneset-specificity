#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 40: Fig 8 (可重复性热图 + 数据集间一致性) + 补充材料表
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

BASE = "."
OUT = os.path.join(BASE, "out")
FIG = os.path.join(BASE, "figures")
SUP = os.path.join(BASE, "supplementary")
os.makedirs(SUP, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "figure.dpi": 300, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.06})

R = pd.read_csv(os.path.join(OUT, "specificity_matrix.csv"))
S = pd.read_csv(os.path.join(OUT, "dataset_summary.csv"))
# 注意：label_permutation.csv 已被 100 次置换的结果覆盖（与
# label_permutation_100perm.csv 逐字节相同），而本脚本需要 200 次置换的列
# （null_median / null_p95 / emp_P）。手稿与 Table S5 用的都是 200 次口径。
LP = pd.read_csv(os.path.join(OUT, "label_permutation_200perm.csv"))
CC = pd.read_csv(os.path.join(OUT, "cell_composition_analysis.csv"))
SURV = json.load(open(os.path.join(OUT, "field_survey_papers.json"), encoding="utf-8"))

DS = ["GSE256068", "GSE213488", "GSE128300"]
wide = R[R.dataset.isin(DS)].pivot_table(index="gene_set", columns="dataset",
                                          values="Z", aggfunc="first")
wide = wide.dropna(thresh=2)
print("[heatmap] gene sets with Z in >=2 datasets:", wide.shape)

# ---------- Fig 8 ----------
fig = plt.figure(figsize=(9.6, 4.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1], wspace=0.34)

# A: heatmap of top variable gene sets
ax = fig.add_subplot(gs[0, 0])
sel = wide.loc[wide.abs().max(axis=1).sort_values(ascending=False).head(40).index]
order = leaves_list(linkage(pdist(sel.fillna(0).values), method="average"))
sel = sel.iloc[order]
vmax = np.nanpercentile(np.abs(sel.values), 98)
im = ax.imshow(sel.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
ax.set_xticks(range(sel.shape[1]))
ax.set_xticklabels([c.replace("GSE", "GSE\n") for c in sel.columns], fontsize=7.5)
ax.set_yticks(range(sel.shape[0]))
ax.set_yticklabels([g[:38] for g in sel.index], fontsize=5.6)
ax.set_title("A  Top-40 most variable gene sets", loc="left",
             fontweight="bold", fontsize=8.8)
cb = fig.colorbar(im, ax=ax, fraction=0.036, pad=0.02)
cb.set_label("Z score", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)

# B: cross-dataset agreement
ax = fig.add_subplot(gs[0, 1])
pairs = [("GSE256068", "GSE213488"), ("GSE256068", "GSE128300"),
         ("GSE213488", "GSE128300")]
cols = ["#2b6cb0", "#c05621", "#2f855a"]
lab = []
for (a, b), c in zip(pairs, cols):
    sub = wide[[a, b]].dropna()
    if len(sub) < 5:
        continue
    r = np.corrcoef(sub[a], sub[b])[0, 1]
    ax.scatter(sub[a], sub[b], s=7, alpha=0.55, color=c,
               label=f"{a.replace('GSE','')} vs {b.replace('GSE','')}  r={r:+.2f} (n={len(sub)})")
    lab.append(r)
ax.axhline(0, color="k", lw=0.6, ls=":"); ax.axvline(0, color="k", lw=0.6, ls=":")
ax.set_xlabel("Gene-set Z score, dataset 1")
ax.set_ylabel("Gene-set Z score, dataset 2")
# 顶部留白：不留白时图例第一行文字被 1 个散点压住（步 134 实测，标记半径 5.5 px
# 时该点圆心正好落在图例文字框内）。
_lo, _hi = ax.get_ylim()
ax.set_ylim(_lo, _hi + 0.18 * (_hi - _lo))
ax.legend(frameon=False, fontsize=6.6, loc="upper left")
ax.set_title("B  Cross-dataset agreement of gene-set Z", loc="left",
             fontweight="bold", fontsize=8.8)
if lab:
    ax.text(0.98, 0.03, "mean r = %+.2f" % np.mean(lab), transform=ax.transAxes,
            ha="right", fontsize=8, fontweight="bold",
            color="#c53030" if np.mean(lab) < 0.3 else "#2f855a")
fig.suptitle("Reproducibility of gene-set enrichment across independent cohorts",
             fontsize=10, y=1.02, x=0.02, ha="left")
fig.savefig(os.path.join(FIG, "Fig8_reproducibility.png"))
plt.close(fig)
print("[Fig8] saved")

# ---------- 补充表 ----------
DS_INFO = {
 "GSE256068": ("Temporal lobe epilepsy with hippocampal sclerosis vs non-epileptic controls",
               "Temporal cortex / hippocampus", "RNA-seq (Illumina NovaSeq 6000)"),
 "GSE213488": ("Focal cortical dysplasia type IIa/IIb vs autopsy controls",
               "Grey matter (cortex)", "RNA-seq (Illumina HiSeq 2500)"),
 "GSE128300": ("Focal cortical dysplasia / HS-epilepsy vs non-epilepsy",
               "Cortex", "RNA-seq (Illumina HiSeq 2500)"),
 "GSE272142": ("Temporal lobe epilepsy vs controls", "Dentate gyrus",
               "RNA-seq (Illumina HiSeq 2500)"),
}
s1 = S.copy()
s1["condition"] = s1.dataset.map(lambda d: DS_INFO.get(d, ("", "", ""))[0])
s1["tissue"] = s1.dataset.map(lambda d: DS_INFO.get(d, ("", "", ""))[1])
s1["platform"] = s1.dataset.map(lambda d: DS_INFO.get(d, ("", "", ""))[2])
lp = LP[["dataset", "null_median", "null_p95", "emp_P"]]
s1 = s1.merge(lp, on="dataset", how="left")
s1.to_csv(os.path.join(SUP, "TableS1_datasets.csv"), index=False)

s2 = R[["dataset", "family", "gene_set", "n", "obs", "exp", "Z", "P", "FDR"]]
s2 = s2.sort_values(["dataset", "P"])
s2.to_csv(os.path.join(SUP, "TableS2_all_genesets.csv"), index=False)

s3 = pd.DataFrame([{"pmid": r["pmid"], "title": r["title"],
                    "reports_specificity_test": r["has_specificity_test"],
                    "methods": ";".join(r["methods"]),
                    "max_reported_DEG": r["max_deg"]} for r in SURV])
s3.to_csv(os.path.join(SUP, "TableS3_field_survey.csv"), index=False)

s4 = CC.copy()
s4.to_csv(os.path.join(SUP, "TableS4_cell_composition.csv"), index=False)

lp.to_csv(os.path.join(SUP, "TableS5_label_permutation.csv"), index=False)

print("\n[supplementary tables]")
for f in sorted(os.listdir(SUP)):
    p = os.path.join(SUP, f)
    n = sum(1 for _ in open(p, encoding="utf-8")) - 1
    print("   %-34s %7.1f KB  %6d rows" % (f, os.path.getsize(p) / 1024, n))

# 关键数字再确认
print("\n[cross-dataset correlation]")
for a, b in pairs:
    sub = wide[[a, b]].dropna()
    if len(sub) >= 5:
        print("   %s vs %s : r = %+.3f (n=%d)" % (a, b, np.corrcoef(sub[a], sub[b])[0, 1], len(sub)))
