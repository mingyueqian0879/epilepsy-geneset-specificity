#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 134a: 为 Fig 8 面板 B 选择图例位置（零遮挡 + 不越界 + 不压标题）。

判定用「图例文字」的真实 bbox（不含 borderpad），并把标记半径算进去。
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

OUT = "out"
R = pd.read_csv(os.path.join(OUT, "specificity_matrix_final.csv"))
wide = R.pivot_table(index="gene_set", columns="dataset", values="Z")
PAIRS = [("GSE256068", "GSE213488"), ("GSE256068", "GSE128300"),
         ("GSE213488", "GSE128300")]
COLS = ["#2b6cb0", "#c05621", "#2f855a"]


def build(figsize, wr, head, fs, ncol, loc, outside=None, csp=1.0):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.linewidth": 0.8, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 300})
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, width_ratios=wr, wspace=0.28)
    ax = fig.add_subplot(gs[0, 1])
    for (a, b), c in zip(PAIRS, COLS):
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
    if head:
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + head * (hi - lo))
    if outside is None:
        leg = ax.legend(frameon=False, fontsize=fs, loc=loc, ncol=ncol,
                        columnspacing=csp)
    else:
        leg = ax.legend(frameon=False, fontsize=fs, loc="upper center",
                        ncol=ncol, bbox_to_anchor=outside, columnspacing=csp)
    ax.set_title("B  Cross-dataset agreement of gene-set Z", loc="left",
                 fontweight="bold", fontsize=8.8)
    ax.text(0.98, 0.03, "mean r = -0.01", transform=ax.transAxes, ha="right",
            fontsize=8, fontweight="bold", color="#c53030")
    fig.canvas.draw()
    return fig, ax, leg


def measure(fig, ax, leg):
    r = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(r) for t in leg.get_texts()]
    pts = []
    for coll in ax.collections:
        offs = coll.get_offsets()
        if offs is None or len(offs) == 0:
            continue
        sizes = coll.get_sizes()
        rad = (np.sqrt(sizes[0]) / 2.0) * (fig.dpi / 72.0) if len(sizes) else 3.0
        pts.append((np.asarray(offs)[:, :2], rad))
    xy = np.vstack([p[0] for p in pts])
    rad = max(p[1] for p in pts)
    disp = ax.transData.transform(xy)
    hit = 0
    for tb in boxes:
        dx = np.maximum(np.maximum(tb.x0 - disp[:, 0], 0), disp[:, 0] - tb.x1)
        dy = np.maximum(np.maximum(tb.y0 - disp[:, 1], 0), disp[:, 1] - tb.y1)
        hit += int((np.hypot(dx, dy) <= rad).sum())
    fb = leg.get_window_extent(r)
    inv = ax.transAxes.inverted()
    (x0, y0), (x1, y1) = inv.transform([[fb.x0, fb.y0], [fb.x1, fb.y1]])
    tb_title = ax.title.get_window_extent(r)
    title_hit = not (fb.y1 < tb_title.y0 or tb_title.y1 < fb.y0
                     or fb.x1 < tb_title.x0 or tb_title.x1 < fb.x0)
    return hit, (x0, y0, x1, y1), rad, title_hit


print("Fig 8 panel B — legend placement scan")
print("  (hit = markers whose ink reaches a legend TEXT box; "
      "fits = legend inside axes horizontally)")
cases = []
for head in (0.0, 0.15, 0.25):
    cases.append(("inside upper-left head=%.2f" % head,
                  dict(head=head, fs=6.6, ncol=1, loc="upper left")))
cases += [
    ("inside lower-left head=0", dict(head=0.0, fs=6.6, ncol=1, loc="lower left")),
    ("inside lower-right head=0", dict(head=0.0, fs=6.6, ncol=1, loc="lower right")),
    ("outside below (1 col)", dict(head=0.0, fs=6.6, ncol=1, loc="upper center",
                                   outside=(0.5, -0.16))),
    ("outside below (2 col)", dict(head=0.0, fs=6.6, ncol=2, loc="upper center",
                                   outside=(0.5, -0.16), csp=1.6)),
    ("inside upper-left head=0.25, fs6.0", dict(head=0.25, fs=6.0, ncol=1,
                                                loc="upper left")),
]
for name, kw in cases:
    try:
        fig, ax, leg = build((9.6, 4.6), [1.35, 1], **kw)
        hit, bb, rad, th = measure(fig, ax, leg)
        fits = bb[0] >= -0.02 and bb[2] <= 1.0
        ok = "OK" if (hit == 0 and fits and not th) else ""
        print("  %-36s hit=%-3d x=[%.2f,%.2f] y=[%.2f,%.2f] fits=%-5s titleClash=%-5s %s"
              % (name, hit, bb[0], bb[2], bb[1], bb[3], fits, th, ok))
        plt.close(fig)
    except Exception as e:
        print("  %-36s ERROR %r" % (name, e))
