#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 113: 生成样本水平 QC 图（补 Fig 2/Fig 3 缺失的 patient-level 证据）
面板：
  A  GSE256068 PCA —— 按分组着色，按组织分形状，直接暴露组×组织混杂
  B  GSE256068 各组织构成（堆叠条）
  C  GSE213488 PCA —— 无混杂的对照
  D  GSE256068 组织分层后的 DEG 率（逐模型）
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "."
M = os.path.join(BASE, "out", "multi")
DST = "figures"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.linewidth": 0.8, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.06})
GRP = {"Case": "#c53030", "Control": "#2b6cb0"}
TIS = {"Temporal": "o", "Frontal": "s", "Parietal": "^", "Hippocampus": "D",
       "grey_matter": "o"}


def pca(X):
    A = X.values.astype(float)
    A = A - A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, keepdims=True); sd[sd == 0] = 1
    A = np.nan_to_num(A / sd)
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    return Vt.T * S, (S ** 2) / (S ** 2).sum()


fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0))

# ---------- A: GSE256068 PCA ----------
X = pd.read_pickle(os.path.join(M, "GSE256068_expr.pkl"))
D = pd.read_csv(os.path.join(M, "GSE256068_design.csv"))
D = D[D["sample"].isin(X.columns)].reset_index(drop=True)
X = X[D["sample"].tolist()]
pc, var = pca(X)
ax = axes[0, 0]
for tis in sorted(D.tissue.unique()):
    for g in ["Case", "Control"]:
        m = (D.tissue == tis) & (D.group == g)
        if not m.any():
            continue
        ax.scatter(pc[m.values, 0], pc[m.values, 1], s=26, marker=TIS.get(tis, "o"),
                   facecolor=GRP[g], edgecolor="k", linewidth=0.4, alpha=0.85,
                   label="%s / %s (n=%d)" % (g, tis, int(m.sum())))
ax.set_xlabel("PC1 (%.1f%%)" % (var[0] * 100))
ax.set_ylabel("PC2 (%.1f%%)" % (var[1] * 100))
ax.axhline(0, color="grey", lw=0.5, ls=":"); ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.legend(fontsize=6.0, frameon=False, loc="upper left", ncol=1, handletextpad=0.3)
ax.set_title("A  GSE256068: group is confounded with brain region",
             loc="left", fontweight="bold", fontsize=8.8)

# ---------- B: 组织构成 ----------
ax = axes[0, 1]
ct = pd.crosstab(D.tissue, D.group)
ct = ct.reindex(columns=[c for c in ["Case", "Control"] if c in ct.columns])
ct.plot(kind="barh", stacked=True, ax=ax,
        color=[GRP[c] for c in ct.columns], width=0.62, edgecolor="k", linewidth=0.4)
for i, tis in enumerate(ct.index):
    tot = ct.loc[tis].sum()
    ax.text(tot + 1.2, i, "n=%d" % int(tot), va="center", fontsize=7)
ax.set_xlabel("Samples")
ax.set_ylabel("")
ax.legend(fontsize=7, frameon=False, loc="lower right")
ax.set_title("B  Only 4 of 27 controls match the cases' region",
             loc="left", fontweight="bold", fontsize=8.8)

# ---------- C: GSE213488 PCA ----------
X2 = pd.read_pickle(os.path.join(M, "GSE213488_expr.pkl"))
D2 = pd.read_csv(os.path.join(M, "GSE213488_design.csv"))
D2 = D2[D2["sample"].isin(X2.columns)].reset_index(drop=True)
X2 = X2[D2["sample"].tolist()]
pc2, var2 = pca(X2)
ax = axes[1, 0]
for g in ["Case", "Control"]:
    m = (D2.group == g).values
    ax.scatter(pc2[m, 0], pc2[m, 1], s=34, marker="o", facecolor=GRP[g],
               edgecolor="k", linewidth=0.4, alpha=0.85, label="%s (n=%d)" % (g, m.sum()))
ax.set_xlabel("PC1 (%.1f%%)" % (var2[0] * 100))
ax.set_ylabel("PC2 (%.1f%%)" % (var2[1] * 100))
ax.axhline(0, color="grey", lw=0.5, ls=":"); ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.legend(fontsize=7, frameon=False)
ax.set_title("C  GSE213488: single tissue, no such confounding",
             loc="left", fontweight="bold", fontsize=8.8)

# ---------- D: 分层 DEG 率 ----------
ax = axes[1, 1]
labels = ["all controls\n(no covariate)", "+ tissue\ncovariate",
          "+ 6 cell scores", "+ tissue\n+ 6 cell scores"]
vals = [59.6, 42.9, 17.0, 0.36]
cols = ["#a0aec0", "#90cdf4", "#63b3ed", "#2b6cb0"]
b = ax.bar(range(4), vals, color=cols, width=0.62, edgecolor="k", linewidth=0.4)
for i, v in enumerate(vals):
    ax.text(i, v + 1.5, "%.1f%%" % v, ha="center", fontsize=7.5)
ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("Genome-wide DEG rate (%)")
ax.set_ylim(0, 70)
ax.set_title("D  GSE256068: adjustment removes almost all signal",
             loc="left", fontweight="bold", fontsize=8.8)

fig.suptitle("Sample-level structure of the cohorts (patient-level evidence absent from Figures 2 and 3)",
             fontsize=10.2, y=1.02, x=0.02, ha="left")
fig.tight_layout()
out = os.path.join(DST, "FigS_sample_structure.png")
fig.savefig(out, dpi=300, facecolor="white")
plt.close(fig)

from PIL import Image
im = Image.open(out)
im.convert("RGB").save(os.path.join(DST, "FigS_sample_structure.tif"), "TIFF",
                       compression="tiff_lzw", dpi=(300, 300))
print("已生成 %s  %dx%d" % (out, im.width, im.height))
