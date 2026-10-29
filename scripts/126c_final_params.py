#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 126c: 锁定最终绘图参数（Fig 2A 标签写法 / Fig 2C 与 Fig S1A 图例布局）"""
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
MATCH = {"GSE256068": "4/27", "GSE213488": "4/4",
         "GSE128300": "8/8", "GSE272142": "6/6"}


def pca2(X):
    A = X.values.astype(float)
    A = A - A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    A = np.nan_to_num(A / sd)
    _, S_, Vt = np.linalg.svd(A, full_matrices=False)
    return (Vt.T * S_)[:, :2], (S_ ** 2) / (S_ ** 2).sum()


def combos_256068():
    X = pd.read_pickle(M + r"\GSE256068_expr.pkl")
    D = pd.read_csv(M + r"\GSE256068_design.csv")
    D = D[D["sample"].isin(X.columns)].reset_index(drop=True)
    pc, var = pca2(X[D["sample"].tolist()])
    out = []
    for t in ["Temporal", "Frontal", "Parietal", "Hippocampus"]:
        for g in ["Case", "Control"]:
            m = ((D.tissue == t) & (D.group == g)).values
            if m.any():
                out.append((g, t, int(m.sum()), m))
    return pc, var, out


print("=== Fig 2A: third-line wording vs clearance (slot ~55-58 px) ===")
S = pd.read_csv("out/dataset_summary.csv").sort_values(
    "deg_rate", ascending=False)
for word in ["region 4/27", "matched 4/27", "match 4/27", "region-matched 4/27"]:
    for fs in (6.0, 6.2):
        labs = ["%s\n(n=%d/%d)\n%s" % (d, c, k, word)
                for d, c, k in zip(S.dataset, S.n_case, S.n_ctrl)]
        fig = plt.figure(figsize=(11.0, 3.4))
        gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15], wspace=0.34)
        ax = fig.add_subplot(gs[0, 0])
        ax.bar(range(len(S)), S.deg_rate,
               color=[PAL[d] for d in S.dataset], width=0.62)
        ax.set_xticks(range(len(S)))
        ax.set_xticklabels(labs, fontsize=fs)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        bbs = [t.get_window_extent(r) for t in ax.get_xticklabels()]
        slot = (ax.transData.transform([(1, 0)])[0][0]
                - ax.transData.transform([(0, 0)])[0][0])
        cl = [(bbs[i + 1].x0 - bbs[i].x1)
              for i in range(len(bbs) - 1)]
        print("  %-22s fs=%.1f slot=%5.1f maxw=%5.1f clearance=%s"
              % (word, fs, slot, max(b.width for b in bbs),
                 " ".join("%+.0f" % c for c in cl)))
        plt.close(fig)

print()
print("=== Fig 2C legend: fit inside axes with 'n=' labels ===")
pc, var, combos = combos_256068()
for fs in (5.6, 5.8, 6.0):
    for csp in (0.4, 0.5, 0.6):
        for head in (0.16, 0.18, 0.20):
            labs = ["%s / %s (n=%d)" % (g, t, n) for g, t, n, _ in combos]
            fig = plt.figure(figsize=(11.0, 3.4))
            gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15], wspace=0.34)
            ax = fig.add_subplot(gs[0, 2])
            for (g, t, n, m), lab in zip(combos, labs):
                ax.scatter(pc[m, 0], pc[m, 1], s=22, marker=TIS[t],
                           facecolor=GRP[g], edgecolor="k", linewidth=0.35,
                           alpha=0.85, label=lab)
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi + head * (hi - lo))
            leg = ax.legend(fontsize=fs, ncol=2, frameon=False, loc="upper left",
                            handletextpad=0.25, columnspacing=csp,
                            labelspacing=0.3, borderpad=0.15)
            fig.canvas.draw()
            r = fig.canvas.get_renderer()
            bb = leg.get_window_extent(r)
            inv = ax.transAxes.inverted()
            (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
            disp = ax.transData.transform(pc)
            m2 = ((disp[:, 0] >= bb.x0) & (disp[:, 0] <= bb.x1)
                  & (disp[:, 1] >= bb.y0) & (disp[:, 1] <= bb.y1))
            ok = "OK" if (m2.sum() == 0 and x1 <= 1.0) else ""
            print("  fs=%.1f csp=%.1f head=%.2f -> occ=%d x1=%.2f %s"
                  % (fs, csp, head, int(m2.sum()), x1, ok))
            plt.close(fig)

print()
print("=== Fig S1A legend: fit inside axes with 'n=' labels ===")
for fs in (6.0, 6.3, 6.5):
    for csp in (0.5, 0.7):
        for head in (0.16, 0.20):
            labs = ["%s / %s (n=%d)" % (g, t, n) for g, t, n, _ in combos]
            fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0))
            ax = axes[0, 0]
            for (g, t, n, m), lab in zip(combos, labs):
                ax.scatter(pc[m, 0], pc[m, 1], s=26, marker=TIS[t],
                           facecolor=GRP[g], edgecolor="k", linewidth=0.4,
                           alpha=0.85, label=lab)
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi + head * (hi - lo))
            leg = ax.legend(fontsize=fs, ncol=2, frameon=False, loc="upper left",
                            handletextpad=0.3, columnspacing=csp, labelspacing=0.3)
            fig.canvas.draw()
            r = fig.canvas.get_renderer()
            bb = leg.get_window_extent(r)
            inv = ax.transAxes.inverted()
            (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
            disp = ax.transData.transform(pc)
            m2 = ((disp[:, 0] >= bb.x0) & (disp[:, 0] <= bb.x1)
                  & (disp[:, 1] >= bb.y0) & (disp[:, 1] <= bb.y1))
            ok = "OK" if (m2.sum() == 0 and x1 <= 1.0) else ""
            print("  fs=%.1f csp=%.1f head=%.2f -> occ=%d x1=%.2f %s"
                  % (fs, csp, head, int(m2.sum()), x1, ok))
            plt.close(fig)
