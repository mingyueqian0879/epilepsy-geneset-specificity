#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 47: 图形摘要 v2 —— 修正文字重叠与溢出"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

PKG = "04_Figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.06, "figure.dpi": 300})

fig, ax = plt.subplots(figsize=(7.6, 4.3), facecolor="white")
ax.set_xlim(0, 100); ax.set_ylim(0, 58); ax.axis("off")

ax.text(1.0, 55.4, "Gene-set associations in epilepsy transcriptomics",
        fontsize=11.5, fontweight="bold", va="center")
ax.text(1.0, 51.2, "do not exceed chance — and the field almost never checks",
        fontsize=10, color="#c53030", fontweight="bold", va="center")

# ---- 三栏 ----
PANELS = [
    (1.0, 31.2, "#2b6cb0", "1  Field practice", "576",
     "epilepsy bioinformatics\npapers (2020–2026)", "only 9 (1.6%) test specificity"),
    (34.9, 30.4, "#c05621", "2  Reanalysis of 4 cohorts", "59.7%",
     "of genes differentially\nexpressed", "0 / 248 sets enriched (FDR)"),
    (68.8, 30.2, "#2f855a", "3  Framework works", "Z = +6.1",
     "microglial markers\nenriched", "real signal IS detected"),
]
for x, w, col, title, big, mid, foot in PANELS:
    ax.add_patch(FancyBboxPatch((x, 20.5), w, 24.5, boxstyle="round,pad=0.5",
                                linewidth=1.2, edgecolor=col, facecolor=col + "12"))
    cx = x + w / 2
    ax.text(cx, 42.4, title, fontsize=8.8, fontweight="bold", color=col,
            ha="center", va="center")
    ax.text(cx, 35.4, big, fontsize=24, fontweight="bold", color=col,
            ha="center", va="center")
    ax.text(cx, 29.6, mid, fontsize=7.3, ha="center", va="center", color="#2d3748")
    ax.text(cx, 23.6, foot, fontsize=7.6, fontweight="bold", color="#c53030",
            ha="center", va="center")

for x0, x1 in [(31.2, 34.6), (65.3, 68.5)]:
    ax.add_patch(FancyArrowPatch((x0, 32.7), (x1, 32.7), arrowstyle="-|>",
                                 mutation_scale=13, color="#a0aec0", lw=1.5))

# ---- 热点基因集 ----
ax.add_patch(FancyBboxPatch((1.0, 11.0), 98, 8.2, boxstyle="round,pad=0.5",
                            linewidth=1.2, edgecolor="#6b46c1", facecolor="#6b46c10e"))
ax.text(50, 16.9, "None of the 8 candidate gene sets reached significance in any cohort",
        fontsize=8.8, fontweight="bold", color="#6b46c1", ha="center", va="center")
ax.text(50, 13.2,
        "lactylation · ferroptosis · cuproptosis · disulfidptosis · pyroptosis · m6A · autophagy · glycolysis",
        fontsize=7.2, ha="center", va="center", color="#2d3748")

# ---- 结论 ----
ax.add_patch(FancyBboxPatch((1.0, 1.2), 98, 8.6, boxstyle="round,pad=0.5",
                            linewidth=1.4, edgecolor="#c53030", facecolor="#c5303010"))
ax.text(50, 7.8, "A gene set being differentially expressed is NOT evidence of involvement",
        fontsize=9.2, fontweight="bold", color="#c53030", ha="center", va="center")
ax.text(50, 4.9, "Signal is real but global, driven by cell-type composition",
        fontsize=7.4, ha="center", va="center", color="#2d3748")
ax.text(50, 2.6, "Cohorts do not agree (mean r = −0.01)  ·  we propose a 5-item reporting standard",
        fontsize=7.4, ha="center", va="center", color="#2d3748")

out = os.path.join(PKG, "graphical_abstract.png")
fig.savefig(out, dpi=300, facecolor="white")
plt.close(fig)
from PIL import Image
im = Image.open(out)
im.save(os.path.join(PKG, "graphical_abstract.tif"), "TIFF",
        compression="tiff_lzw", dpi=(300, 300))
print("[graphical abstract v2] %dx%d px  PNG %.0f KB  TIFF %.0f KB" % (
    im.width, im.height, os.path.getsize(out) / 1024,
    os.path.getsize(os.path.join(PKG, "graphical_abstract.tif")) / 1024))
