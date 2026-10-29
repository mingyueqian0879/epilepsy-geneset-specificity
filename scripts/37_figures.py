#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 37: 生成论文配图 (publication-ready, English labels)
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = "."
OUT = os.path.join(BASE, "out")
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621",
       "GSE128300": "#2f855a", "GSE272142": "#6b46c1"}

R = pd.read_csv(os.path.join(OUT, "specificity_matrix.csv"))
S = pd.read_csv(os.path.join(OUT, "dataset_summary.csv"))
SURV = json.load(open(os.path.join(OUT, "field_survey_papers.json"), encoding="utf-8"))

# ============ Fig 1: 领域调查 ============
fig, axes = plt.subplots(1, 3, figsize=(9.5, 2.9))
n = len(SURV)
nsp = sum(1 for r in SURV if r["has_specificity_test"])
n_deg = sum(1 for r in SURV if r["max_deg"])

ax = axes[0]
ax.bar([0, 1], [n - nsp, nsp], color=["#cbd5e0", "#e53e3e"], width=0.6)
ax.set_xticks([0, 1]); ax.set_xticklabels(["No specificity\ntesting", "Reports\nspecificity test"])
ax.set_ylabel("Published papers (n)")
for x, v in zip([0, 1], [n - nsp, nsp]):
    ax.text(x, v + 8, f"{v}\n({100*v/n:.1f}%)", ha="center", va="bottom", fontsize=8)
ax.set_ylim(0, n * 1.15); ax.set_title("A  Field practice", loc="left", fontweight="bold")

ax = axes[1]
from collections import Counter
c = Counter()
for r in SURV:
    for m in r["methods"]:
        c[m] += 1
top = c.most_common(7)[::-1]
ax.barh([t[0] for t in top], [100 * t[1] / n for t in top], color="#4299e1")
ax.set_xlabel("% of papers"); ax.set_title("B  Methods used", loc="left", fontweight="bold")
for i, t in enumerate(top):
    ax.text(100 * t[1] / n + 0.1, i, f"{t[1]}", va="center", fontsize=7)

ax = axes[2]
degs = [r["max_deg"] for r in SURV if r["max_deg"]]
ax.hist(np.log10(degs), bins=18, color="#805ad5", edgecolor="white", linewidth=0.4)
ax.set_xlabel("log10(reported DEG count)"); ax.set_ylabel("Papers")
ax.axvline(np.log10(np.median(degs)), color="k", ls="--", lw=1)
ax.text(np.log10(np.median(degs)), ax.get_ylim()[1] * 0.9,
        f" median={int(np.median(degs))}", fontsize=7)
ax.set_title("C  Reported DEG counts", loc="left", fontweight="bold")
fig.suptitle(f"Survey of {n} epilepsy bioinformatics papers (2020-2026)",
             fontsize=9.5, y=1.06, x=0.02, ha="left")
fig.savefig(os.path.join(FIG, "Fig1_field_survey.png"))
plt.close(fig)
print("[Fig1] saved")

# ============ Fig 2: DEG 率 ============
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.1))
ax = axes[0]
ds = S.sort_values("deg_rate", ascending=False)
cols = [PAL.get(d, "#888") for d in ds.dataset]
b = ax.bar(range(len(ds)), ds.deg_rate, color=cols, width=0.62)
ax.set_xticks(range(len(ds)))
ax.set_xticklabels([f"{d}\n(n={c}/{k})" for d, c, k in
                    zip(ds.dataset, ds.n_case, ds.n_ctrl)], fontsize=7.5)
ax.set_ylabel("DEG rate at FDR<0.05 (%)")
ax.axhspan(5, 20, color="gray", alpha=0.18, zorder=0)
ax.text(len(ds) - 0.5, 12.5, "typical\n5-20%", fontsize=7, ha="right", color="#4a5568")
for i, v in enumerate(ds.deg_rate):
    ax.text(i, v + 1.2, f"{v:.1f}%", ha="center", fontsize=8)
ax.set_ylim(0, 70)
ax.set_title("A  Same clinical question, DEG rate spans 0-60%",
             loc="left", fontweight="bold", fontsize=8.5)

ax = axes[1]
labels, vals, errs = [], [], []
for d in ds.dataset:
    sub = R[R.dataset == d]
    if not len(sub):
        continue
    labels.append(d); vals.append(sub.Z.median())
    errs.append([sub.Z.median() - sub.Z.quantile(.25), sub.Z.quantile(.75) - sub.Z.median()])
if labels:
    ax.bar(range(len(labels)), vals, yerr=np.array(errs).T, capsize=3,
           color=[PAL.get(l, "#888") for l in labels], width=0.6)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([l.replace("GSE", "GSE\n") for l in labels], fontsize=7.5)
    ax.set_ylabel("Gene-set Z score (median, IQR)")
    ax.set_title("B  Gene-set enrichment centered at zero",
                 loc="left", fontweight="bold", fontsize=8.5)
fig.savefig(os.path.join(FIG, "Fig2_deg_and_Z.png"))
plt.close(fig)
print("[Fig2] saved")

# ============ Fig 3: Z 分布小提琴 ============
fig, ax = plt.subplots(figsize=(5.6, 3.2))
data, labs, cols2 = [], [], []
for d in ["GSE256068", "GSE213488", "GSE128300"]:
    sub = R[R.dataset == d].dropna(subset=["Z"])
    if len(sub) > 5:
        data.append(sub.Z.values); labs.append(f"{d}\n(n={len(sub)})"); cols2.append(PAL[d])
if data:
    vp = ax.violinplot(data, showmedians=True, widths=0.75)
    for i, b in enumerate(vp["bodies"]):
        b.set_facecolor(cols2[i]); b.set_alpha(0.55); b.set_edgecolor("k"); b.set_linewidth(0.6)
    for k in ("cmedians", "cbars", "cmins", "cmaxes"):
        if k in vp: vp[k].set_color("k"); vp[k].set_linewidth(0.8)
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axhline(2, color="#e53e3e", lw=0.8, ls=":")
    ax.axhline(-2, color="#e53e3e", lw=0.8, ls=":")
    ax.text(0.62, 2.1, "|Z|=2", color="#e53e3e", fontsize=7)
    ax.set_xticks(range(1, len(labs) + 1)); ax.set_xticklabels(labs, fontsize=7.5)
    ax.set_ylabel("Gene-set enrichment Z score")
    ax.set_title("Distribution of Z scores across 248 gene sets\n(MSigDB Hallmark + KEGG)",
                 loc="left", fontweight="bold", fontsize=8.5)
fig.savefig(os.path.join(FIG, "Fig3_Z_distribution.png"))
plt.close(fig)
print("[Fig3] saved")

# ============ Fig 4: 热点基因集 vs 阳性对照 ============
fig, ax = plt.subplots(figsize=(7.4, 3.6))
hot = R[R.family == "HOT"].copy()
order = (hot.groupby("gene_set").Z.max().sort_values(ascending=True).index.tolist())
ds_order = [d for d in ["GSE256068", "GSE213488"] if d in set(hot.dataset)]
w = 0.36
for j, d in enumerate(ds_order):
    vals, pos = [], []
    for i, g in enumerate(order):
        sub = hot[(hot.gene_set == g) & (hot.dataset == d)]
        vals.append(sub.Z.iloc[0] if len(sub) else np.nan)
        pos.append(i + (j - 0.5) * w)
    ax.barh(pos, vals, height=w, color=PAL[d], label=d, alpha=0.9)
ax.axvline(0, color="k", lw=0.8)
ax.axvline(1.96, color="#e53e3e", lw=0.9, ls=":")
ax.text(2.05, -0.8, "Z=1.96", color="#e53e3e", fontsize=7)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([g.replace("_core", "").replace("_regulators", "").replace("_LRG332", "")
                    for g in order], fontsize=8)
ax.set_xlabel("Gene-set enrichment Z score")
ax.legend(frameon=False, fontsize=7.5, loc="lower right")
ax.set_title("Candidate 'hot' gene sets vs cell-type positive controls",
             loc="left", fontweight="bold", fontsize=9)
for i, g in enumerate(order):
    if g.startswith("POS_"):
        ax.get_yticklabels()[i].set_fontweight("bold")
fig.savefig(os.path.join(FIG, "Fig4_hot_vs_positive.png"))
plt.close(fig)
print("[Fig4] saved")

print("\n[figures in %s]" % FIG)
for f in sorted(os.listdir(FIG)):
    print("   %-34s %6.0f KB" % (f, os.path.getsize(os.path.join(FIG, f)) / 1024))
