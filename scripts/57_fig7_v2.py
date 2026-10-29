#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 57: 重画 Figure 7 —— 报告标准由 5 条改为 6 条"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

FIG = "04_Figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.bbox": "tight",
                     "savefig.pad_inches": 0.06, "figure.dpi": 300})

STEPS = [
    ("1  Differential expression",
     "Report the GLOBAL DEG rate\n(not just hits inside the gene set)", "#2b6cb0"),
    ("2  Specificity test   REQUIRED",
     "Expression-matched permutation\n\u2192 report Z and empirical P", "#c05621"),
    ("3  Multiple-testing correction",
     "Report FDR-adjusted results\n(not nominal P alone)", "#2f855a"),
    ("4  Detectable effect size   NEW",
     "State the smallest enrichment the design could detect\n(a null result is otherwise uninterpretable)", "#805ad5"),
    ("5  Positive control",
     "Report as a calibration probe (rank, effect size),\nand state whether it entered the FDR family", "#6b46c1"),
    ("6  Composition adjustment",
     "Report DEGs before AND after\ncell-composition adjustment", "#b7791f"),
]

fig, ax = plt.subplots(figsize=(7.9, 5.1))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

ax.text(1.0, 96.0, "Minimum reporting standard for gene-set studies",
        fontsize=11, fontweight="bold", va="center")

y = 88.5
box_h = 12.6
gap = 1.4
for title, body, col in STEPS:
    ax.add_patch(FancyBboxPatch((1.0, y - box_h), 98, box_h,
                                boxstyle="round,pad=0.4",
                                linewidth=1.15, edgecolor=col, facecolor=col + "18"))
    ax.text(3.0, y - 3.4, title, fontsize=9.0, fontweight="bold", color=col, va="center")
    ax.text(3.0, y - 8.4, body, fontsize=7.4, color="#2d3748", va="center", linespacing=1.45)
    y -= (box_h + gap)

for i in range(len(STEPS) - 1):
    y0 = 88.5 - i * (box_h + gap) - box_h
    ax.add_patch(FancyArrowPatch((50, y0 - 0.15), (50, y0 - gap + 0.15),
                                 arrowstyle="-|>", mutation_scale=10,
                                 color="#a0aec0", lw=1.0))

ax.text(1.0, 1.6, "A gene set being differentially expressed is NOT evidence of involvement.",
        fontsize=8.6, style="italic", color="#c53030", va="center")

out = os.path.join(FIG, "Fig7_reporting_standard.png")
fig.savefig(out, dpi=300, facecolor="white")
plt.close(fig)
from PIL import Image
im = Image.open(out)
im.convert("RGB").save(os.path.join(FIG, "Fig7_reporting_standard.tif"), "TIFF",
                       compression="tiff_lzw", dpi=(300, 300))
print("[Fig7 v2] %dx%d px   PNG %.0f KB   TIFF %.0f KB" % (
    im.width, im.height, os.path.getsize(out) / 1024,
    os.path.getsize(os.path.join(FIG, "Fig7_reporting_standard.tif")) / 1024))
