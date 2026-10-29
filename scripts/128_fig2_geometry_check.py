#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 128: 对重绘后的 Fig 2 做「标题 vs 图例」「图例 vs 散点」的精确几何复核，
并把图例区域裁成放大图供目视确认。

判定全部在 fig.dpi 对应的显示坐标系里做，单位是像素，不受 bbox_inches='tight' 影响。
"""
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

M = "out/multi"
OUT = "out/_fig2_verify"
os.makedirs(OUT, exist_ok=True)
GRP = {"Case": "#c53030", "Control": "#2b6cb0"}
TIS = {"Temporal": "o", "Frontal": "s", "Parietal": "^", "Hippocampus": "D"}
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621",
       "GSE128300": "#2f855a", "GSE272142": "#6b46c1"}


def rects_intersect(a, b):
    return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])


fig = plt.figure(figsize=(11.0, 3.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15], wspace=0.34)

# A
S = pd.read_csv("out/dataset_summary.csv").sort_values(
    "deg_rate", ascending=False)
MATCH = {"GSE256068": "4/27", "GSE213488": "4/4",
         "GSE128300": "8/8", "GSE272142": "6/6"}
axA = fig.add_subplot(gs[0, 0])
axA.bar(range(len(S)), S.deg_rate, color=[PAL[d] for d in S.dataset], width=0.62)
axA.set_xticks(range(len(S)))
axA.set_xticklabels(["%s\n(n=%d/%d)\nregion %s" % (d, c, k, MATCH[d])
                     for d, c, k in zip(S.dataset, S.n_case, S.n_ctrl)], fontsize=6.0)
axA.set_ylabel("DEG rate at FDR<0.05 (%)")
axA.set_ylim(0, 70)
axA.set_title("A  Same clinical question, DEG rate spans 0-60%", loc="left",
              fontweight="bold", fontsize=8.2)

# B (geometry irrelevant, kept so panel C sits in the same place)
R = pd.read_csv("out/specificity_matrix_final.csv")
axB = fig.add_subplot(gs[0, 1])
labs, vals, errs = [], [], []
for d in S.dataset:
    sub = R[R.dataset == d].dropna(subset=["Z"])
    if not len(sub):
        continue
    labs.append(d); vals.append(sub.Z.median())
    errs.append([sub.Z.median() - sub.Z.quantile(.25),
                 sub.Z.quantile(.75) - sub.Z.median()])
axB.bar(range(len(labs)), vals, yerr=np.array(errs).T, capsize=3,
        color=[PAL[l] for l in labs], width=0.6)
axB.set_xticks(range(len(labs)))
axB.set_xticklabels([l.replace("GSE", "GSE\n") for l in labs], fontsize=6.6)
axB.set_ylabel("Gene-set Z score (median, IQR)")
axB.set_title("B  Enrichment centred at zero", loc="left", fontweight="bold",
              fontsize=8.2)

# C
X = pd.read_pickle(os.path.join(M, "GSE256068_expr.pkl"))
D = pd.read_csv(os.path.join(M, "GSE256068_design.csv"))
D = D[D["sample"].isin(X.columns)].reset_index(drop=True)
X = X[D["sample"].tolist()]
A = X.values.astype(float)
A = A - A.mean(axis=1, keepdims=True)
sd = A.std(axis=1, keepdims=True); sd[sd == 0] = 1
A = np.nan_to_num(A / sd)
_, S_, Vt = np.linalg.svd(A, full_matrices=False)
pc = (Vt.T * S_)[:, :2]
var = (S_ ** 2) / (S_ ** 2).sum()

axC = fig.add_subplot(gs[0, 2])
for tis in ["Temporal", "Frontal", "Parietal", "Hippocampus"]:
    for g in ["Case", "Control"]:
        m = ((D.tissue == tis) & (D.group == g)).values
        if not m.any():
            continue
        axC.scatter(pc[m, 0], pc[m, 1], s=22, marker=TIS[tis], facecolor=GRP[g],
                    edgecolor="k", linewidth=0.35, alpha=0.85,
                    label="%s / %s (n=%d)" % (g, tis, int(m.sum())))
axC.set_xlabel("PC1 (%.1f%%)" % (var[0] * 100), fontsize=7.5)
axC.set_ylabel("PC2 (%.1f%%)" % (var[1] * 100), fontsize=7.5)
axC.tick_params(labelsize=6.5)
lo, hi = axC.get_ylim()
axC.set_ylim(lo, hi + 0.18 * (hi - lo))
leg = axC.legend(fontsize=5.8, ncol=2, frameon=False, loc="upper left",
                 handletextpad=0.25, columnspacing=0.4, labelspacing=0.3,
                 borderpad=0.15)
axC.set_title("C  GSE256068 patient-level PCA:\ngroup is confounded with brain region",
              loc="left", fontweight="bold", fontsize=8.2)

fig.canvas.draw()
r = fig.canvas.get_renderer()

lb = leg.get_window_extent(r)
tb = axC.title.get_window_extent(r)
ab = axC.get_window_extent(r)

print("legend  bbox (px): x=[%.1f,%.1f] y=[%.1f,%.1f]" % (lb.x0, lb.x1, lb.y0, lb.y1))
print("title   bbox (px): x=[%.1f,%.1f] y=[%.1f,%.1f]" % (tb.x0, tb.x1, tb.y0, tb.y1))
print("axes    bbox (px): x=[%.1f,%.1f] y=[%.1f,%.1f]" % (ab.x0, ab.x1, ab.y0, ab.y1))

print("\n[overlap checks]")
print("  legend vs title  : %s"
      % ("OVERLAP" if rects_intersect((lb.x0, lb.y0, lb.x1, lb.y1),
                                      (tb.x0, tb.y0, tb.x1, tb.y1)) else "clear"))
print("  legend vs axes   : %s"
      % ("OVERLAP" if rects_intersect((lb.x0, lb.y0, lb.x1, lb.y1),
                                      (ab.x0, ab.y0, ab.x1, ab.y1)) else "clear"))

# strict point test in PIXEL space, with marker radius
radius = (np.sqrt(22) / 2.0) * (fig.dpi / 72.0)
disp = axC.transData.transform(pc[:, :2])
d = np.sqrt(
    np.maximum(np.maximum(lb.x0 - disp[:, 0], disp[:, 0] - lb.x1), 0) ** 2 +
    np.maximum(np.maximum(lb.y0 - disp[:, 1], disp[:, 1] - lb.y1), 0) ** 2)
n_hit = int((d <= radius).sum())
print("  points within marker radius (%.1f px) of legend rect: %d" % (radius, n_hit))
print("  closest point centre distance to legend rect: %.1f px" % d.min())

# zoomed crop of the legend region (with margin)
pad = 40
x0 = max(int(lb.x0) - pad, 0); x1 = min(int(lb.x1) + pad, int(fig.bbox.x1))
y0 = max(int(lb.y0) - pad, 0); y1 = min(int(lb.y1) + pad, int(fig.bbox.y1))
fig.savefig(os.path.join(OUT, "fig2_full.png"), dpi=300, facecolor="white")
plt.close(fig)

im = Image.open(os.path.join(OUT, "fig2_full.png"))
h = im.height
# display coords have origin at bottom; image origin is top
crop = im.crop((x0, h - y1, x1, h - y0))
crop = crop.resize((crop.width * 3, crop.height * 3), Image.LANCZOS)
crop.save(os.path.join(OUT, "fig2C_legend_zoom.png"))
print("\nzoomed legend crop -> %s  (%dx%d)"
      % (os.path.join(OUT, "fig2C_legend_zoom.png"), crop.width, crop.height))
