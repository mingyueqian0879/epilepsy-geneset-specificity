#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 126b: 两个版面问题的量化诊断
(1) Fig 2A 的 x 轴三行标签是否在水平方向互相重叠（渲染宽度 vs 刻度间距）
(2) Fig 2C / Fig S1A 图例在若干标签写法下的实际 bbox，选出「既零遮挡又不越界」的方案
"""
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = "out/multi"
GRP = {"Case": "#c53030", "Control": "#2b6cb0"}
TIS = {"Temporal": "o", "Frontal": "s", "Parietal": "^", "Hippocampus": "D"}
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621",
       "GSE128300": "#2f855a", "GSE272142": "#6b46c1"}


def pca(X):
    A = X.values.astype(float)
    A = A - A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    A = np.nan_to_num(A / sd)
    _, S_, Vt = np.linalg.svd(A, full_matrices=False)
    return (Vt.T * S_)[:, :2], (S_ ** 2) / (S_ ** 2).sum()


def load(ds):
    X = pd.read_pickle(os.path.join(M, ds + "_expr.pkl"))
    D = pd.read_csv(os.path.join(M, ds + "_design.csv"))
    D = D[D["sample"].isin(X.columns)].reset_index(drop=True)
    return X[D["sample"].tolist()], D


def sample_cols(D):
    out = []
    for tis in ["Temporal", "Frontal", "Parietal", "Hippocampus"]:
        for g in ["Case", "Control"]:
            m = ((D.tissue == tis) & (D.group == g)).values
            if m.any():
                out.append((g, tis, int(m.sum()), m))
    return out


# ---------------- (1) Fig 2A tick-label overlap ----------------
def check_2A():
    print("=" * 74)
    print("Fig 2A x tick labels: rendered width vs slot width")
    print("=" * 74)
    S = pd.read_csv("out/dataset_summary.csv")
    S = S.sort_values("deg_rate", ascending=False)
    MATCH = {"GSE256068": "4/27", "GSE213488": "4/4",
             "GSE128300": "8/8", "GSE272142": "6/6"}
    fig = plt.figure(figsize=(11.0, 3.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15], wspace=0.34)
    ax = fig.add_subplot(gs[0, 0])
    ax.bar(range(len(S)), S.deg_rate,
           color=[PAL.get(d, "#888") for d in S.dataset], width=0.62)
    ax.set_xticks(range(len(S)))
    labs = ["%s\n(n=%d/%d)\nregion matched: %s" %
            (d, c, k, MATCH.get(d, "?"))
            for d, c, k in zip(S.dataset, S.n_case, S.n_ctrl)]
    ax.set_xticklabels(labs, fontsize=6.4)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    print("  slot width (1 data unit) = %.1f px" % (
        ax.transData.transform([(1, 0)])[0][0] - ax.transData.transform([(0, 0)])[0][0]))
    for i, t in enumerate(ax.get_xticklabels()):
        bb = t.get_window_extent(r)
        print("  [%d] width=%6.1f px  text=%r" % (i, bb.width, labs[i].replace("\n", " | ")))
    # horizontal overlap between adjacent labels
    bbs = [t.get_window_extent(r) for t in ax.get_xticklabels()]
    for i in range(len(bbs) - 1):
        ov = bbs[i].x1 - bbs[i + 1].x0
        print("  overlap %d-%d : %+6.1f px %s" % (
            i, i + 1, ov, "<== OVERLAP" if ov > 0 else ""))
    plt.close(fig)
    print()


# ---------------- (2) legend variants ----------------
def legend_bbox(ax, pts, labels, fs, ncol, head, loc="upper left"):
    if head:
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + head * (hi - lo))
    handles = ax.get_legend_handles_labels()[0]
    leg = ax.legend(handles, labels, fontsize=fs, frameon=False, loc=loc, ncol=ncol,
                    handletextpad=0.3, labelspacing=0.3, borderpad=0.15,
                    columnspacing=0.9)
    fig = ax.figure
    fig.canvas.draw()
    bb = leg.get_window_extent(fig.canvas.get_renderer())
    inv = ax.transAxes.inverted()
    (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    disp = ax.transData.transform(np.asarray(pts)[:, :2])
    m = ((disp[:, 0] >= bb.x0) & (disp[:, 0] <= bb.x1) &
         (disp[:, 1] >= bb.y0) & (disp[:, 1] <= bb.y1))
    fits = (x0 >= -0.02) and (x1 <= 1.0)
    return int(m.sum()), (x0, y0, x1, y1), fits


def variant_test():
    print("=" * 74)
    print("Legend variants (goal: occluded=0 AND bbox fits inside axes)")
    print("=" * 74)
    X, D = load("GSE256068")
    pc, var = pca(X)
    combos = sample_cols(D)
    variants = {
        "long":   ["%s / %s (n=%d)" % (g, t, n) for g, t, n, _ in combos],
        "nodot":  ["%s / %s (%d)" % (g, t, n) for g, t, n, _ in combos],
    }
    for fs in (5.8, 6.0, 6.4):
        for ncol in (2,):
            for head in (0.18, 0.22, 0.26):
                for vname, labs in variants.items():
                    fig = plt.figure(figsize=(11.0, 3.4))
                    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15],
                                          wspace=0.34)
                    ax = fig.add_subplot(gs[0, 2])
                    for (g, t, n, m), lab in zip(combos, labs):
                        ax.scatter(pc[m, 0], pc[m, 1], s=22, marker=TIS[t],
                                   facecolor=GRP[g], edgecolor="k", linewidth=0.35,
                                   alpha=0.85, label=lab)
                    ax.set_xlabel("PC1 (%.1f%%)" % (var[0] * 100), fontsize=7.5)
                    ax.set_ylabel("PC2 (%.1f%%)" % (var[1] * 100), fontsize=7.5)
                    n_occ, bb, fits = legend_bbox(ax, pc, labs, fs, ncol, head)
                    ok = "OK" if (n_occ == 0 and fits) else ""
                    print("  %-7s fs=%.1f ncol=%d head=%.2f -> occ=%d "
                          "bbox=(%.2f,%.2f,%.2f,%.2f) fits=%s %s"
                          % (vname, fs, ncol, head, n_occ, bb[0], bb[1], bb[2],
                             bb[3], fits, ok))
                    plt.close(fig)
    print()


if __name__ == "__main__":
    check_2A()
    variant_test()
