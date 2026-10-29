#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 126: 图例遮挡的确定性诊断（精确 bbox 计算，不依赖目测）。

复现 Fig 2C / Fig S1A 的当前绘图代码，计算「图例矩形」覆盖了多少个散点；
再遍历候选方案（位置 x 列数 x 字号 x y 轴余量），输出每个方案的遮挡数。
遮挡数 = 0 才是安全方案。
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


def pca(X):
    A = X.values.astype(float)
    A = A - A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    A = np.nan_to_num(A / sd)
    _, S_, Vt = np.linalg.svd(A, full_matrices=False)
    return Vt.T * S_, (S_ ** 2) / (S_ ** 2).sum()


def load(ds):
    X = pd.read_pickle(os.path.join(M, ds + "_expr.pkl"))
    D = pd.read_csv(os.path.join(M, ds + "_design.csv"))
    D = D[D["sample"].isin(X.columns)].reset_index(drop=True)
    X = X[D["sample"].tolist()]
    return X, D


def occluded(ax, leg, pts):
    """返回 (被图例矩形盖住的散点数, 图例在 axes 坐标下的 bbox)。"""
    fig = ax.figure
    fig.canvas.draw()
    bb = leg.get_window_extent(fig.canvas.get_renderer())
    inv = ax.transAxes.inverted()
    (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    disp = ax.transData.transform(np.asarray(pts)[:, :2])
    m = ((disp[:, 0] >= bb.x0) & (disp[:, 0] <= bb.x1) &
         (disp[:, 1] >= bb.y0) & (disp[:, 1] <= bb.y1))
    return int(m.sum()), (x0, y0, x1, y1)


def scatter_256068(ax, pt_size=22):
    X, D = load("GSE256068")
    pc, var = pca(X)
    n_lab = 0
    for tis in ["Temporal", "Frontal", "Parietal", "Hippocampus"]:
        for g in ["Case", "Control"]:
            m = ((D.tissue == tis) & (D.group == g)).values
            if not m.any():
                continue
            n_lab += 1
            ax.scatter(pc[m, 0], pc[m, 1], s=pt_size, marker=TIS.get(tis, "o"),
                       facecolor=GRP[g], edgecolor="k", linewidth=0.35,
                       alpha=0.85,
                       label="%s / %s (n=%d)" % (g, tis, int(m.sum())))
    ax.set_xlabel("PC1 (%.1f%%)" % (var[0] * 100), fontsize=7.5)
    ax.set_ylabel("PC2 (%.1f%%)" % (var[1] * 100), fontsize=7.5)
    ax.axhline(0, color="grey", lw=0.5, ls=":")
    ax.axvline(0, color="grey", lw=0.5, ls=":")
    ax.tick_params(labelsize=6.5)
    return pc, n_lab


def run(tag, make_fig, cases, pt_size):
    print("=" * 74)
    print(tag)
    print("=" * 74)
    best = None
    for name, loc, ncol, fs, head in cases:
        fig, ax = make_fig()
        pts, n_lab = scatter_256068(ax, pt_size=pt_size)
        if head:
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi + head * (hi - lo))
        leg = ax.legend(fontsize=fs, frameon=False, loc=loc, ncol=ncol,
                        handletextpad=0.3, labelspacing=0.3, borderpad=0.15,
                        columnspacing=1.0)
        n, bb = occluded(ax, leg, pts)
        flag = "  <== SAFE" if n == 0 else ""
        print("  %-42s occluded=%2d  bbox=(%.2f,%.2f,%.2f,%.2f)%s"
              % (name, n, bb[0], bb[1], bb[2], bb[3], flag))
        if best is None or n < best[1]:
            best = (name, n)
        plt.close(fig)
    print("  legend entries rendered: %d" % n_lab)
    print("  best: %s  (occluded=%d)" % best)
    print()
    return best


def fig2c():
    def make():
        fig = plt.figure(figsize=(11.0, 3.4))
        gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15], wspace=0.34)
        return fig, fig.add_subplot(gs[0, 2])
    cases = [
        ("current: upper left, 1col, fs5.2", "upper left", 1, 5.2, 0.0),
        ("lower left, 1col, fs5.2", "lower left", 1, 5.2, 0.0),
        ("upper left, 2col, fs5.2", "upper left", 2, 5.2, 0.0),
        ("lower left, 2col, fs5.2", "lower left", 2, 5.2, 0.0),
        ("lower left, 2col, fs6.0", "lower left", 2, 6.0, 0.0),
        ("upper left, 2col, fs6.0, +20% y", "upper left", 2, 6.0, 0.20),
        ("upper left, 2col, fs6.0, +35% y", "upper left", 2, 6.0, 0.35),
        ("upper left, 1col, fs6.0, +45% y", "upper left", 1, 6.0, 0.45),
    ]
    return run("Fig 2C  (11.0 x 3.4 in; panel C = GSE256068 PCA)", make, cases, 22)


def figs1a():
    def make():
        fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0))
        return fig, axes[0, 0]
    cases = [
        ("current: upper left, 1col, fs6.0", "upper left", 1, 6.0, 0.0),
        ("lower left, 1col, fs6.0", "lower left", 1, 6.0, 0.0),
        ("upper left, 2col, fs6.0", "upper left", 2, 6.0, 0.0),
        ("lower left, 2col, fs6.0", "lower left", 2, 6.0, 0.0),
        ("lower left, 2col, fs6.5", "lower left", 2, 6.5, 0.0),
        ("upper left, 2col, fs6.5, +20% y", "upper left", 2, 6.5, 0.20),
        ("upper left, 2col, fs6.5, +32% y", "upper left", 2, 6.5, 0.32),
        ("upper left, 1col, fs6.5, +42% y", "upper left", 1, 6.5, 0.42),
    ]
    return run("Fig S1A  (9.6 x 7.0 in 2x2; panel A = GSE256068 PCA)", make, cases, 26)


if __name__ == "__main__":
    a = fig2c()
    b = figs1a()
    print("SUMMARY")
    print("  Fig 2C best : %s -> occluded=%d" % a)
    print("  Fig S1A best: %s -> occluded=%d" % b)
