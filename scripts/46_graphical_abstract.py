#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 46: 图形摘要 (graphical abstract) —— BIB 要求: 单独文件、landscape、命名 graphical_abstract
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

PKG = "04_Figures"
os.makedirs(PKG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.06, "figure.dpi": 300})

fig, ax = plt.subplots(figsize=(7.4, 4.15), facecolor="white")
ax.set_xlim(0, 100); ax.set_ylim(0, 56); ax.axis("off")

# ---- 标题 ----
ax.text(1.5, 52.6, "Gene-set associations in epilepsy transcriptomics",
        fontsize=11.5, fontweight="bold", va="center")
ax.text(1.5, 48.6, "do not exceed chance — and the field almost never checks",
        fontsize=10, color="#c53030", fontweight="bold", va="center")

# ---- 左: 领域现状 ----
box = FancyBboxPatch((1.5, 27.5), 29, 17.5, boxstyle="round,pad=0.6",
                     linewidth=1.2, edgecolor="#2b6cb0", facecolor="#2b6cb012")
ax.add_patch(box)
ax.text(16, 42.4, "1  Field practice", fontsize=9.5, fontweight="bold",
        color="#2b6cb0", ha="center", va="center")
ax.text(16, 36.4, "576", fontsize=26, fontweight="bold", color="#2b6cb0",
        ha="center", va="center")
ax.text(16, 31.9, "epilepsy bioinformatics papers\n(2020–2026)", fontsize=7.6,
        ha="center", va="center", color="#2d3748")
ax.text(16, 29.4, "only 9 (1.6%) test specificity", fontsize=8.2, fontweight="bold",
        color="#c53030", ha="center", va="center")

# ---- 中: 核心结果 ----
box2 = FancyBboxPatch((34, 27.5), 32, 17.5, boxstyle="round,pad=0.6",
                      linewidth=1.2, edgecolor="#c05621", facecolor="#c0562112")
ax.add_patch(box2)
ax.text(50, 42.4, "2  Reanalysis of 4 cohorts", fontsize=9.5, fontweight="bold",
        color="#c05621", ha="center", va="center")
ax.text(50, 36.6, "59.7%", fontsize=26, fontweight="bold", color="#c05621",
        ha="center", va="center")
ax.text(50, 32.2, "of genes differentially expressed", fontsize=7.6,
        ha="center", va="center", color="#2d3748")
ax.text(50, 29.4, "yet 0 / 248 gene sets enriched (FDR)", fontsize=8.2,
        fontweight="bold", color="#c53030", ha="center", va="center")

# ---- 右: 方法有效 ----
box3 = FancyBboxPatch((69.5, 27.5), 29, 17.5, boxstyle="round,pad=0.6",
                      linewidth=1.2, edgecolor="#2f855a", facecolor="#2f855a12")
ax.add_patch(box3)
ax.text(84, 42.4, "3  Framework works", fontsize=9.5, fontweight="bold",
        color="#2f855a", ha="center", va="center")
ax.text(84, 36.6, "Z = +6.1", fontsize=26, fontweight="bold", color="#2f855a",
        ha="center", va="center")
ax.text(84, 32.2, "microglial markers enriched", fontsize=7.6,
        ha="center", va="center", color="#2d3748")
ax.text(84, 29.4, "real signal IS detected", fontsize=8.2, fontweight="bold",
        color="#2f855a", ha="center", va="center")

# ---- 箭头 ----
for x0, x1 in [(30.6, 33.8), (66.2, 69.4)]:
    ax.add_patch(FancyArrowPatch((x0, 36.2), (x1, 36.2), arrowstyle="-|>",
                                 mutation_scale=14, color="#a0aec0", lw=1.6))

# ---- 底部: 热点基因集全灭 ----
box4 = FancyBboxPatch((1.5, 15.5), 97, 9.2, boxstyle="round,pad=0.6",
                      linewidth=1.2, edgecolor="#6b46c1", facecolor="#6b46c10e")
ax.add_patch(box4)
ax.text(50, 22.3, "None of the 8 candidate gene sets reached significance in any cohort",
        fontsize=9, fontweight="bold", color="#6b46c1", ha="center", va="center")
ax.text(50, 18.3,
        "lactylation   ·   ferroptosis   ·   cuproptosis   ·   disulfidptosis   ·   "
        "pyroptosis   ·   m6A   ·   autophagy   ·   glycolysis",
        fontsize=8, ha="center", va="center", color="#2d3748")

# ---- 结论条 ----
box5 = FancyBboxPatch((1.5, 4.2), 97, 9.5, boxstyle="round,pad=0.6",
                      linewidth=1.4, edgecolor="#c53030", facecolor="#c5303010")
ax.add_patch(box5)
ax.text(50, 10.6, "A gene set being differentially expressed is NOT evidence of involvement",
        fontsize=9.6, fontweight="bold", color="#c53030", ha="center", va="center")
ax.text(50, 6.6,
        "The signal is real but global, driven by cell-type composition  ·  "
        "enrichment scores uncorrelated between cohorts (mean r = −0.01)\n"
        "We propose a 5-item minimum reporting standard and release the framework",
        fontsize=7.7, ha="center", va="center", color="#2d3748")

fig.savefig(os.path.join(PKG, "graphical_abstract.png"), dpi=300, facecolor="white")
from PIL import Image
im = Image.open(os.path.join(PKG, "graphical_abstract.png"))
im.save(os.path.join(PKG, "graphical_abstract.tif"), "TIFF",
        compression="tiff_lzw", dpi=(300, 300))
print("[graphical abstract] %dx%d px   PNG %.0f KB   TIFF %.0f KB" % (
    im.width, im.height,
    os.path.getsize(os.path.join(PKG, "graphical_abstract.png")) / 1024,
    os.path.getsize(os.path.join(PKG, "graphical_abstract.tif")) / 1024))
