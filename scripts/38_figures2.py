#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 38: 剩余配图 (Fig 5-7)"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = "."
OUT = os.path.join(BASE, "out")
FIG = os.path.join(BASE, "figures")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.05})
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621",
       "GSE128300": "#2f855a", "GSE272142": "#6b46c1"}

# 注意：label_permutation.csv 已被 100 次置换的结果覆盖；本脚本的 Fig 5 需要
# 200 次置换的列（null_median / null_p95 / emp_P）。
LP = pd.read_csv(os.path.join(OUT, "label_permutation_200perm.csv"))
CC = pd.read_csv(os.path.join(OUT, "cell_composition_analysis.csv"))

# ============ Fig 5: 标签置换校准 ============
fig, axes = plt.subplots(1, 4, figsize=(11, 2.9), sharey=False)
for ax, (_, r) in zip(axes, LP.iterrows()):
    ds = r.dataset
    f = os.path.join(OUT, f"labelperm_{ds}.rds")
    # RDS 不可直接读 -> 用已存的汇总值绘制示意
    ax.bar([0], [r.obs_deg], color=PAL.get(ds, "#888"), width=0.55, label="Observed")
    ax.bar([1], [r.null_median], color="#a0aec0", width=0.55, label="Null (median)")
    ax.bar([1], [r.null_p95], color="#e2e8f0", width=0.55, label="Null 95th")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Observed", "Label-\npermuted"], fontsize=7.5)
    ax.set_title(f"{ds}\n(n={r.n_case}/{r.n_ctrl})", fontsize=8, fontweight="bold")
    ax.text(0, r.obs_deg * 1.03 + 2, f"{int(r.obs_deg):,}", ha="center", fontsize=7.5)
    ax.text(1, r.null_p95 * 1.03 + 2, f"p95={r.null_p95:.0f}", ha="center", fontsize=7)
    lab = ("p = %.3f" % r.emp_P) + ("\n(real signal)" if r.emp_P < 0.05 else "\n(NO signal)")
    ax.text(0.5, 0.80, lab, transform=ax.transAxes, ha="center", fontsize=8,
            color="#c53030" if r.emp_P >= 0.05 else "#2f855a", fontweight="bold")
    ax.set_ylabel("DEG count (FDR<0.05)" if ax is axes[0] else "")
fig.suptitle("Label-permutation calibration: is the differential-expression signal real?",
             fontsize=9.5, y=1.10, x=0.02, ha="left")
fig.savefig(os.path.join(FIG, "Fig5_label_permutation.png"))
plt.close(fig)
print("[Fig5] saved")

# ============ Fig 6: 细胞组成校正 ============
deg = CC[CC.marker.str.startswith("DEG_")].copy()
piv = deg.pivot_table(index="dataset", columns="marker", values="obs", aggfunc="first")
piv = piv.reindex(["GSE256068", "GSE213488", "GSE128300", "GSE272142"]).dropna(how="all")

fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.1))
ax = axes[0]
x = np.arange(len(piv)); w = 0.36
ax.bar(x - w/2, piv["DEG_modelA"], w, label="Model A: group only", color="#e53e3e", alpha=0.85)
ax.bar(x + w/2, piv["DEG_modelB_celladj"], w, label="Model B: + cell-composition", color="#2b6cb0", alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(piv.index, fontsize=7.5)
ax.set_ylabel("DEG count (FDR<0.05)")
ax.legend(frameon=False, fontsize=7.5)
for i, (a, b) in enumerate(zip(piv["DEG_modelA"], piv["DEG_modelB_celladj"])):
    if a > 50:
        ax.text(i, a * 1.04, f"-{100*(a-b)/a:.0f}%", ha="center", fontsize=7.5,
                color="#c53030", fontweight="bold")
ax.set_title("A  Cell-composition adjustment\nremoves most DEGs", loc="left",
             fontweight="bold", fontsize=8.5)

ax = axes[1]
mk = CC[~CC.marker.str.startswith("DEG_")].copy()
mk["is_pos"] = mk.marker.isin(["micro", "astro", "oligo", "neuro", "endo", "immune"])
p = mk.pivot_table(index="marker", columns="dataset", values="Z", aggfunc="first")
p = p.dropna(how="all")
if len(p):
    yy = np.arange(len(p))
    for j, d in enumerate([c for c in p.columns if c in PAL]):
        ax.barh(yy + (j - 0.5) * 0.35, p[d].values, height=0.35,
                color=PAL[d], label=d)
    ax.set_yticks(yy); ax.set_yticklabels(p.index, fontsize=8)
    ax.axvline(0, color="k", lw=0.8)
    ax.legend(frameon=False, fontsize=7)
ax.set_xlabel("Cell-marker gene-set Z score")
ax.set_title("B  Cell-type markers are strongly\nenriched (pipeline works)",
             loc="left", fontweight="bold", fontsize=8.5)
fig.savefig(os.path.join(FIG, "Fig6_cell_composition.png"))
plt.close(fig)
print("[Fig6] saved")

# ============ Fig 7: 建议流程 ============
fig, ax = plt.subplots(figsize=(7.6, 4.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
steps = [
    (8.6, "1  Differential expression", "Report the GLOBAL DEG rate\n(not just hits inside the gene set)", "#2b6cb0"),
    (6.9, "2  Specificity test  REQUIRED", "Expression-matched permutation\n-> report Z and empirical P", "#c05621"),
    (5.2, "3  Multiple-testing correction", "Report FDR-adjusted results\n(not nominal P alone)", "#2f855a"),
    (3.5, "4  Positive control", "Include a cell-type marker set\nto prove the pipeline can detect real signal", "#6b46c1"),
    (1.8, "5  Composition adjustment", "Report DEGs before AND after\ncell-composition adjustment", "#b7791f"),
]
for y, title, body, col in steps:
    box = FancyBboxPatch((0.4, y - 0.62), 9.2, 1.28, boxstyle="round,pad=0.08",
                         linewidth=1.1, edgecolor=col, facecolor=col + "18")
    ax.add_patch(box)
    ax.text(0.75, y + 0.30, title, fontsize=9.5, fontweight="bold", color=col, va="center")
    ax.text(0.75, y - 0.22, body, fontsize=8, color="#2d3748", va="center")
for y0, y1 in zip([8.6, 6.9, 5.2, 3.5], [6.9, 5.2, 3.5, 1.8]):
    ax.add_patch(FancyArrowPatch((5.0, y0 - 0.68), (5.0, y1 + 0.72),
                                 arrowstyle="-|>", mutation_scale=12, color="#a0aec0", lw=1.1))
ax.text(0.4, 9.55, "Minimum reporting standard for epilepsy gene-set studies",
        fontsize=10.5, fontweight="bold")
ax.text(0.4, 0.55, "A gene set being differentially expressed is NOT evidence of involvement.",
        fontsize=9, style="italic", color="#c53030")
fig.savefig(os.path.join(FIG, "Fig7_reporting_standard.png"))
plt.close(fig)
print("[Fig7] saved")

print("\n[figures]")
for f in sorted(os.listdir(FIG)):
    print("   %-36s %6.0f KB" % (f, os.path.getsize(os.path.join(FIG, f)) / 1024))
