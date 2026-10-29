#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 133: Fig 8 面板 B 的「图例 vs 散点」精确判定。

上一版按图例整体 bbox（含 borderpad 内边距）判定，可能高估。
这里改为逐个「图例文字」的真实 bbox 判定，并给出像素级间隙与落盘坐标，
以判断该接触在印刷尺度上是否可见。
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = "."
OUT = os.path.join(BASE, "out")
TMP = os.path.join(OUT, "_layout_repro")
SUBFIG = "figures/Fig8_reproducibility.png"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.linewidth": 0.8, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.06})

R = pd.read_csv(os.path.join(OUT, "specificity_matrix_final.csv"))
wide = R.pivot_table(index="gene_set", columns="dataset", values="Z")

fig = plt.figure(figsize=(9.6, 4.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], wspace=0.28)
ax = fig.add_subplot(gs[0, 1])

pairs = [("GSE256068", "GSE213488"), ("GSE256068", "GSE128300"),
         ("GSE213488", "GSE128300")]
cols = ["#2b6cb0", "#c05621", "#2f855a"]
for (a, b), c in zip(pairs, cols):
    sub = wide[[a, b]].dropna()
    if len(sub) < 5:
        continue
    r = np.corrcoef(sub[a], sub[b])[0, 1]
    ax.scatter(sub[a], sub[b], s=7, alpha=0.55, color=c,
               label="%s vs %s  r=%+.2f (n=%d)"
                     % (a.replace("GSE", ""), b.replace("GSE", ""), r, len(sub)))
ax.axhline(0, color="k", lw=0.6, ls=":")
ax.axvline(0, color="k", lw=0.6, ls=":")
ax.set_xlabel("Gene-set Z score, dataset 1")
ax.set_ylabel("Gene-set Z score, dataset 2")
leg = ax.legend(frameon=False, fontsize=6.6, loc="upper left")
ax.set_title("B  Cross-dataset agreement of gene-set Z", loc="left",
             fontweight="bold", fontsize=8.8)
ax.text(0.98, 0.03, "mean r = -0.01", transform=ax.transAxes, ha="right",
        fontsize=8, fontweight="bold", color="#c53030")

fig.canvas.draw()
r_ = fig.canvas.get_renderer()
dpi_scale = 300.0 / fig.dpi  # display px -> saved-image px

print("figure dpi = %.0f ; saved image is %.2fx display px" % (fig.dpi, dpi_scale))

# all scatter points in display coords
pts = []
for coll in ax.collections:
    offs = coll.get_offsets()
    if offs is None or len(offs) == 0:
        continue
    sizes = coll.get_sizes()
    rad = (np.sqrt(sizes[0]) / 2.0) * (fig.dpi / 72.0) if len(sizes) else 3.0
    pts.append((np.asarray(offs)[:, :2], rad))
pts_xy = np.vstack([p[0] for p in pts])
rad = max(p[1] for p in pts)
print("scatter points: %d ; marker radius = %.1f px (display) = %.1f px (saved)"
      % (len(pts_xy), rad, rad * dpi_scale))

disp = ax.transData.transform(pts_xy)

leg_box = leg.get_window_extent(r_)
print("\nlegend frame bbox (display px): x=[%.1f,%.1f] y=[%.1f,%.1f]"
      % (leg_box.x0, leg_box.x1, leg_box.y0, leg_box.y1))
print("  note: frame includes borderpad; text boxes below are the visible ink\n")

worst = None
for i, t in enumerate(leg.get_texts()):
    tb = t.get_window_extent(r_)
    d = np.sqrt(
        np.maximum(np.maximum(tb.x0 - disp[:, 0], 0), disp[:, 0] - tb.x1) ** 2 +
        np.maximum(np.maximum(tb.y0 - disp[:, 1], 0), disp[:, 1] - tb.y1) ** 2)
    k = int(np.argmin(d))
    touch = int((d <= rad).sum())
    print("  legend text %d %-34r  text bbox x=[%.0f,%.0f] y=[%.0f,%.0f]"
          % (i, t.get_text()[:34], tb.x0, tb.x1, tb.y0, tb.y1))
    print("      nearest marker: data=(%.2f, %.2f)  gap = %+.1f px ; markers within radius: %d"
          % (pts_xy[k, 0], pts_xy[k, 1], d[k], touch))
    if worst is None or d[k] < worst[0]:
        worst = (d[k], i, pts_xy[k], tb)

print("\nWORST case: legend line %d, gap = %+.1f px (display) = %+.1f px (at 300 dpi)"
      % (worst[1], worst[0], worst[0] * dpi_scale))
print("  marker radius at 300 dpi = %.1f px" % (rad * dpi_scale))
verdict = ("VISIBLE overlap (marker ink reaches the glyphs)"
           if worst[0] < rad else
           "no ink-level overlap: nearest glyph box is clear of every marker")
print("  verdict: %s" % verdict)

n_inset = int(((disp[:, 0] >= leg_box.x0) & (disp[:, 0] <= leg_box.x1)
               & (disp[:, 1] >= leg_box.y0) & (disp[:, 1] <= leg_box.y1)).sum())
print("  markers whose CENTRE lies inside the legend frame: %d" % n_inset)

# where is it in the saved image? map display -> saved px
fig.savefig(os.path.join(TMP, "_fig8_panelB_probe.png"), dpi=300, facecolor="white")
scale = 300.0 / fig.dpi
print("\nfor cropping the delivered PNG, the legend occupies approx "
      "x=[%.0f,%.0f] y=[%.0f,%.0f] in the figure canvas (display px)"
      % (leg_box.x0, leg_box.x1, leg_box.y0, leg_box.y1))
print("  -> multiply by %.2f and measure from the TOP for image coordinates" % scale)
