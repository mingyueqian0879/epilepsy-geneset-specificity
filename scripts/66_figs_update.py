#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 66: 用修正后的 P/FDR 重画 Fig 3、Fig 4，并新增 Fig 9 (spiked-in 验证曲线)"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "."
OUT = os.path.join(BASE, "out")
FIG = "04_Figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.06})
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621", "GSE128300": "#2f855a"}

R = pd.read_csv(os.path.join(OUT, "specificity_matrix_final.csv"))
SP = pd.read_csv(os.path.join(OUT, "spikein_simulation_corrected.csv"))

# ---------- Fig 3 (重画: Z 分布, 用最终数据) ----------
fig, ax = plt.subplots(figsize=(5.8, 3.3))
data, labs, cols = [], [], []
for d in ["GSE256068", "GSE213488"]:
    s = R[R.dataset == d].dropna(subset=["Z"])
    if len(s) > 5:
        data.append(s.Z.values); labs.append("%s\n(n=%d)" % (d, len(s))); cols.append(PAL[d])
vp = ax.violinplot(data, showmedians=True, widths=0.72)
for i, b in enumerate(vp["bodies"]):
    b.set_facecolor(cols[i]); b.set_alpha(0.55); b.set_edgecolor("k"); b.set_linewidth(0.6)
for k in ("cmedians", "cbars", "cmins", "cmaxes"):
    if k in vp: vp[k].set_color("k"); vp[k].set_linewidth(0.8)
ax.axhline(0, color="k", lw=0.8, ls="--")
ax.axhline(2, color="#e53e3e", lw=0.8, ls=":")
ax.axhline(-2, color="#e53e3e", lw=0.8, ls=":")
ax.text(1.55, 2.15, "|Z| = 2", color="#e53e3e", fontsize=7)
ax.set_xticks(range(1, len(labs) + 1)); ax.set_xticklabels(labs, fontsize=7.5)
ax.set_ylabel("Gene-set enrichment Z score")
ax.set_title("Distribution of Z scores across 248 gene sets\n(MSigDB Hallmark + KEGG)",
             loc="left", fontweight="bold", fontsize=8.6)
fig.savefig(os.path.join(FIG, "Fig3_Z_distribution.png"), dpi=300, facecolor="white")
plt.close(fig)
print("[Fig3] re-saved")

# ---------- Fig 4 (重画: 热点集 vs 校准探针, 加 FDR 线) ----------
fig, ax = plt.subplots(figsize=(7.6, 3.8))
hot = R[R.family == "HOT"].copy()
order = (hot.groupby("gene_set").Z.max().sort_values(ascending=True).index.tolist())
ds_order = [d for d in ["GSE256068", "GSE213488"] if d in set(hot.dataset)]
w = 0.36
for j, d in enumerate(ds_order):
    vals, pos = [], []
    for i, g in enumerate(order):
        s = hot[(hot.gene_set == g) & (hot.dataset == d)]
        vals.append(s.Z.iloc[0] if len(s) else np.nan)
        pos.append(i + (j - 0.5) * w)
    ax.barh(pos, vals, height=w, color=PAL[d], label=d, alpha=0.9)
ax.axvline(0, color="k", lw=0.8)
ax.axvline(1.96, color="#e53e3e", lw=0.9, ls=":")
ax.text(2.05, -0.85, "Z = 1.96", color="#e53e3e", fontsize=7)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([g.replace("_core", "").replace("_regulators", "")
                    .replace("_LRG332", "").replace("POS_", "") for g in order], fontsize=8)
ax.set_xlabel("Gene-set enrichment Z score")
ax.legend(frameon=False, fontsize=7.5, loc="lower right")
ax.set_title("Candidate 'hot' gene sets vs cell-type calibration probes",
             loc="left", fontweight="bold", fontsize=9)
for i, g in enumerate(order):
    if g.startswith("POS_"):
        ax.get_yticklabels()[i].set_fontweight("bold")
fig.savefig(os.path.join(FIG, "Fig4_hot_vs_positive.png"), dpi=300, facecolor="white")
plt.close(fig)
print("[Fig4] re-saved")

# ---------- Fig 9 (新增: spiked-in 验证) ----------
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2))

ax = axes[0]
ax.plot(SP.median_excess_pct, SP.detection_rate, "o-", color="#2b6cb0", lw=1.6, ms=5)
ax.axhline(90, color="#a0aec0", ls=":", lw=0.9)
ax.axvline(31, color="#c53030", ls="--", lw=1.1)
ax.text(31.6, 45, "analytic MDE\n(+31%)", color="#c53030", fontsize=7)
ax.set_xlabel("Achieved excess of DEGs in spiked set (%)")
ax.set_ylabel("Detection rate (%)")
ax.set_ylim(-5, 105); ax.set_xlim(0, 50)
ax.set_title("A  Detection rate rises only\nfor large enrichment", loc="left",
             fontweight="bold", fontsize=8.6)
for _, r in SP.iterrows():
    if r.detection_rate > 0 or r.median_excess_pct == 1.5:
        ax.annotate("%.0f%%" % r.detection_rate, (r.median_excess_pct, r.detection_rate),
                    textcoords="offset points", xytext=(3, 6), fontsize=7)

ax = axes[1]
ax.plot(SP.median_excess_pct, SP.median_Z, "s-", color="#c05621", lw=1.6, ms=5, label="observed")
ax.axhline(3.72, color="#c53030", ls="--", lw=1.1, label="BH threshold (Z = 3.72)")
ax.set_xlabel("Achieved excess of DEGs in spiked set (%)")
ax.set_ylabel("Median Z score")
ax.set_xlim(0, 50)
ax.legend(frameon=False, fontsize=7.5, loc="upper left")
ax.set_title("B  Null calibration: Z = 0.18 at δ = 0\n(false-positive rate 0%)",
             loc="left", fontweight="bold", fontsize=8.6)
ax.text(42, 0.6, "null", fontsize=7, color="#4a5568", ha="center")

fig.suptitle("Spiked-in simulation validates the expression-matched null",
             fontsize=9.8, y=1.05, x=0.02, ha="left")
fig.savefig(os.path.join(FIG, "Fig9_spikein_validation.png"), dpi=300, facecolor="white")
plt.close(fig)
print("[Fig9] created")

# ---------- TIFF 副本 ----------
from PIL import Image
for f in ["Fig3_Z_distribution", "Fig4_hot_vs_positive", "Fig9_spikein_validation"]:
    im = Image.open(os.path.join(FIG, f + ".png"))
    im.convert("RGB").save(os.path.join(FIG, f + ".tif"), "TIFF",
                           compression="tiff_lzw", dpi=(300, 300))
    print("   %-32s %dx%d  PNG %.0fKB  TIFF %.0fKB" % (
        f, im.width, im.height,
        os.path.getsize(os.path.join(FIG, f + ".png")) / 1024,
        os.path.getsize(os.path.join(FIG, f + ".tif")) / 1024))
