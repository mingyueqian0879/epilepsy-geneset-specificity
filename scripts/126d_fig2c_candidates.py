#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 126d: Fig 2C 图例的两种候选布局对比（生成对照图，供人工/视觉复核）
候选 1：图例留在面板内，2 列 + y 轴留白
候选 2：图例移到面板外（正下方，2 列），与 A 面板三行标签的深度对齐
"""
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = "out/multi"
OUT = "out/_fig2c_candidates"
os.makedirs(OUT, exist_ok=True)
GRP = {"Case": "#c53030", "Control": "#2b6cb0"}
TIS = {"Temporal": "o", "Frontal": "s", "Parietal": "^", "Hippocampus": "D"}
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621",
       "GSE128300": "#2f855a", "GSE272142": "#6b46c1"}
MATCH = {"GSE256068": "4/27", "GSE213488": "4/4",
         "GSE128300": "8/8", "GSE272142": "6/6"}


def pca2(X):
    A = X.values.astype(float)
    A = A - A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    A = np.nan_to_num(A / sd)
    _, S_, Vt = np.linalg.svd(A, full_matrices=False)
    return (Vt.T * S_)[:2].T, (S_ ** 2) / (S_ ** 2).sum()


def base_panels(fig):
    S = pd.read_csv("out/dataset_summary.csv")
    S = S.sort_values("deg_rate", ascending=False)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15], wspace=0.34)

    # A
    ax = fig.add_subplot(gs[0, 0])
    ax.bar(range(len(S)), S.deg_rate,
           color=[PAL.get(d, "#888") for d in S.dataset], width=0.62)
    ax.set_xticks(range(len(S)))
    ax.set_xticklabels(["%s\n(n=%d/%d)\nregion %s"
                        % (d, c, k, MATCH.get(d, "?"))
                        for d, c, k in zip(S.dataset, S.n_case, S.n_ctrl)],
                       fontsize=6.0)
    ax.set_ylabel("DEG rate at FDR<0.05 (%)")
    ax.axhspan(5, 20, color="gray", alpha=0.18, zorder=0)
    ax.text(len(S) - 0.5, 12.5, "typical\n5-20%", fontsize=7, ha="right",
            color="#4a5568")
    for i, v in enumerate(S.deg_rate):
        ax.text(i, v + 1.2, "%.1f%%" % v, ha="center", fontsize=8)
    ax.set_ylim(0, 70)
    ax.set_title("A  Same clinical question, DEG rate spans 0-60%", loc="left",
                 fontweight="bold", fontsize=8.2)

    # B
    R = pd.read_csv("out/specificity_matrix_final.csv")
    ax = fig.add_subplot(gs[0, 1])
    labels, vals, errs = [], [], []
    for d in S.dataset:
        sub = R[R.dataset == d].dropna(subset=["Z"])
        if not len(sub):
            continue
        labels.append(d)
        vals.append(sub.Z.median())
        errs.append([sub.Z.median() - sub.Z.quantile(.25),
                     sub.Z.quantile(.75) - sub.Z.median()])
    ax.bar(range(len(labels)), vals, yerr=np.array(errs).T, capsize=3,
           color=[PAL.get(l, "#888") for l in labels], width=0.6)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([l.replace("GSE", "GSE\n") for l in labels], fontsize=6.6)
    ax.set_ylabel("Gene-set Z score (median, IQR)")
    ax.set_title("B  Enrichment centred at zero", loc="left", fontweight="bold",
                 fontsize=8.2)

    # C
    X = pd.read_pickle(os.path.join(M, "GSE256068_expr.pkl"))
    D = pd.read_csv(os.path.join(M, "GSE256068_design.csv"))
    D = D[D["sample"].isin(X.columns)].reset_index(drop=True)
    X = X[D["sample"].tolist()]
    pc, var = pca2(X)
    ax = fig.add_subplot(gs[0, 2])
    for tis in ["Temporal", "Frontal", "Parietal", "Hippocampus"]:
        for g in ["Case", "Control"]:
            m = ((D.tissue == tis) & (D.group == g)).values
            if not m.any():
                continue
            ax.scatter(pc[m, 0], pc[m, 1], s=22, marker=TIS.get(tis, "o"),
                       facecolor=GRP[g], edgecolor="k", linewidth=0.35,
                       alpha=0.85,
                       label="%s / %s (n=%d)" % (g, tis, int(m.sum())))
    ax.set_xlabel("PC1 (%.1f%%)" % (var[0] * 100), fontsize=7.5)
    ax.set_ylabel("PC2 (%.1f%%)" % (var[1] * 100), fontsize=7.5)
    ax.axhline(0, color="grey", lw=0.5, ls=":")
    ax.axvline(0, color="grey", lw=0.5, ls=":")
    ax.tick_params(labelsize=6.5)
    ax.set_title("C  GSE256068 patient-level PCA:\ngroup is confounded with brain region",
                 loc="left", fontweight="bold", fontsize=8.2)
    return ax, pc


def cand1():
    """Legend inside, 2 columns, y headroom."""
    fig = plt.figure(figsize=(11.0, 3.4))
    ax, pc = base_panels(fig)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.18 * (hi - lo))
    ax.legend(fontsize=5.8, ncol=2, frameon=False, loc="upper left",
              handletextpad=0.25, columnspacing=0.4, labelspacing=0.3,
              borderpad=0.15)
    fig.canvas.draw()
    leg = ax.get_legend()
    bb = leg.get_window_extent(fig.canvas.get_renderer())
    inv = ax.transAxes.inverted()
    (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    disp = ax.transData.transform(pc)
    m = ((disp[:, 0] >= bb.x0) & (disp[:, 0] <= bb.x1)
         & (disp[:, 1] >= bb.y0) & (disp[:, 1] <= bb.y1))
    print("  cand1 legend-inside : occluded=%d  x=[%.2f,%.2f] y=[%.2f,%.2f]"
          % (int(m.sum()), x0, x1, y0, y1))
    fig.savefig(os.path.join(OUT, "cand1_inside.png"), dpi=300, facecolor="white")
    plt.close(fig)


def cand2():
    """Legend outside, directly below panel C, 2 columns."""
    fig = plt.figure(figsize=(11.0, 3.4))
    ax, pc = base_panels(fig)
    leg = ax.legend(fontsize=6.2, ncol=2, frameon=False, loc="upper center",
                    bbox_to_anchor=(0.5, -0.30), handletextpad=0.3,
                    columnspacing=1.4, labelspacing=0.4)
    fig.canvas.draw()
    bb = leg.get_window_extent(fig.canvas.get_renderer())
    inv = ax.transAxes.inverted()
    (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    print("  cand2 legend-below  : bbox x=[%.2f,%.2f] y=[%.2f,%.2f] (outside axes, "
          "occlusion impossible)" % (x0, x1, y0, y1))
    fig.savefig(os.path.join(OUT, "cand2_below.png"), dpi=300, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    cand1()
    cand2()
    print("  written to", OUT)
