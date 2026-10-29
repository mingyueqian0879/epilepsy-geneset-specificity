#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 39: 重做 Fig 5 —— 标签置换校准 (修正布局)"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "."
OUT = os.path.join(BASE, "out")
FIG = os.path.join(BASE, "figures")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.06})
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621",
       "GSE128300": "#2f855a", "GSE272142": "#6b46c1"}

# 注意：label_permutation.csv 已被 100 次置换的结果覆盖；本脚本需要 200 次
# 置换的列（null_median / null_p95 / emp_P）。手稿 Figure 5B 亦为 200 次口径。
LP = pd.read_csv(os.path.join(OUT, "label_permutation_200perm.csv"))
LP["obs_pct"] = 100 * LP.obs_deg / LP.n_genes
LP["null_pct"] = 100 * LP.null_median / LP.n_genes
LP["null_p95_pct"] = 100 * LP.null_p95 / LP.n_genes

fig = plt.figure(figsize=(9.4, 3.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], wspace=0.30)

# ---- Panel A: DEG rate, observed vs null (log scale) ----
ax = fig.add_subplot(gs[0, 0])
x = np.arange(len(LP)); w = 0.34
obs = LP.obs_pct.clip(lower=0.01)
nul = LP.null_p95_pct.clip(lower=0.01)
ax.bar(x - w/2, obs, w, color=[PAL[d] for d in LP.dataset], label="Observed", zorder=3)
ax.bar(x + w/2, nul, w, color="#cbd5e0", label="Null (95th pct)", zorder=3)
ax.set_yscale("log")
ax.set_ylim(0.008, 300)
ax.set_xticks(x)
ax.set_xticklabels([f"{d}\n({c}/{k})" for d, c, k in
                    zip(LP.dataset, LP.n_case, LP.n_ctrl)], fontsize=7.5)
ax.set_ylabel("DEG rate, % (log scale)")
for i, (o, n) in enumerate(zip(obs, nul)):
    ax.text(i - w/2, o * 1.5, f"{LP.obs_pct.iloc[i]:.1f}%", ha="center", fontsize=7.5,
            fontweight="bold", color=PAL[LP.dataset.iloc[i]])
    ax.text(i + w/2, n * 1.5, f"{LP.null_p95_pct.iloc[i]:.2f}%", ha="center", fontsize=7,
            color="#4a5568")
ax.legend(frameon=False, fontsize=7.5, loc="upper right", ncol=1)
ax.axhline(5, color="k", ls=":", lw=0.8)
ax.text(len(LP) - 0.42, 5.6, "5%", fontsize=7, color="#4a5568")
ax.set_title("A  Observed DEG rate vs label-permuted null",
             loc="left", fontweight="bold", fontsize=8.8)
ax.grid(axis="y", alpha=0.25, zorder=0)

# ---- Panel B: empirical P ----
ax = fig.add_subplot(gs[0, 1])
cols = ["#2f855a" if p < 0.05 else "#c53030" for p in LP.emp_P]
b = ax.barh(np.arange(len(LP)), LP.emp_P, color=cols, height=0.6, zorder=3)
ax.set_yticks(np.arange(len(LP)))
ax.set_yticklabels(LP.dataset, fontsize=8)
ax.axvline(0.05, color="k", ls="--", lw=1)
ax.text(0.052, -0.62, "P = 0.05", fontsize=7.5)
for i, (p, d, o) in enumerate(zip(LP.emp_P, LP.dataset, LP.obs_deg)):
    tag = "real signal" if p < 0.05 else "NO signal"
    ax.text(min(p, 1.0) + 0.02, i, f"P = {p:.3f}  ({tag})",
            va="center", fontsize=7.5,
            color="#2f855a" if p < 0.05 else "#c53030", fontweight="bold")
ax.set_xlim(0, 1.45)
ax.set_xlabel("Empirical P (observed DEG vs 200 label permutations)")
ax.set_title("B  Is there any differential-expression signal?",
             loc="left", fontweight="bold", fontsize=8.8)
ax.grid(axis="x", alpha=0.25, zorder=0)

fig.suptitle("Label-permutation calibration of the differential-expression step",
             fontsize=10, y=1.03, x=0.02, ha="left")
fig.savefig(os.path.join(FIG, "Fig5_label_permutation.png"))
plt.close(fig)
print("[Fig5] re-saved")
print("   大小 = %.0f KB" % (os.path.getsize(os.path.join(FIG, "Fig5_label_permutation.png")) / 1024))
