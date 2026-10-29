#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 114: 重绘 Fig 2（新增 C 面板 = 样本水平 PCA）+ 生成 Supplementary Figure S1
P0-3：把队列混杂画进主图，使其不可能被忽略
P2  ：Fig 2A 的 x 轴标注各队列的脑区匹配度
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "."
OUT = os.path.join(BASE, "out")
M = os.path.join(OUT, "multi")
DST = "figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.linewidth": 0.8, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.06})
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621",
       "GSE128300": "#2f855a", "GSE272142": "#6b46c1"}
GRP = {"Case": "#c53030", "Control": "#2b6cb0"}
TIS = {"Temporal": "o", "Frontal": "s", "Parietal": "^", "Hippocampus": "D",
       "grey_matter": "o", "cortex": "o", "dentate_gyrus": "o"}

R = pd.read_csv(os.path.join(OUT, "specificity_matrix_final.csv"))
S = pd.read_csv(os.path.join(OUT, "dataset_summary.csv"))

# 各队列的脑区匹配度（对照中有多少例与病例同组织）
MATCH = {"GSE256068": "4/27", "GSE213488": "4/4", "GSE128300": "8/8", "GSE272142": "6/6"}


def pca(X):
    A = X.values.astype(float)
    A = A - A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, keepdims=True); sd[sd == 0] = 1
    A = np.nan_to_num(A / sd)
    _, S_, Vt = np.linalg.svd(A, full_matrices=False)
    return Vt.T * S_, (S_ ** 2) / (S_ ** 2).sum()


# ================= Fig 2（三面板） =================
fig = plt.figure(figsize=(11.0, 3.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15], wspace=0.34)

# --- A: DEG 率 + 脑区匹配度 ---
ax = fig.add_subplot(gs[0, 0])
ds = S.sort_values("deg_rate", ascending=False)
ax.bar(range(len(ds)), ds.deg_rate, color=[PAL.get(d, "#888") for d in ds.dataset], width=0.62)
ax.set_xticks(range(len(ds)))
ax.set_xticklabels(["%s\n(n=%d/%d)\nregion matched: %s" %
                    (d, c, k, MATCH.get(d, "?"))
                    for d, c, k in zip(ds.dataset, ds.n_case, ds.n_ctrl)], fontsize=6.4)
ax.set_ylabel("DEG rate at FDR<0.05 (%)")
ax.axhspan(5, 20, color="gray", alpha=0.18, zorder=0)
ax.text(len(ds) - 0.5, 12.5, "typical\n5-20%", fontsize=7, ha="right", color="#4a5568")
for i, v in enumerate(ds.deg_rate):
    ax.text(i, v + 1.2, "%.1f%%" % v, ha="center", fontsize=8)
ax.set_ylim(0, 70)
ax.set_title("A  Same clinical question, DEG rate spans 0-60%", loc="left",
             fontweight="bold", fontsize=8.2)

# --- B: 中位 Z ---
ax = fig.add_subplot(gs[0, 1])
labels, vals, errs = [], [], []
for d in ds.dataset:
    sub = R[R.dataset == d].dropna(subset=["Z"])
    if not len(sub):
        continue
    labels.append(d); vals.append(sub.Z.median())
    errs.append([sub.Z.median() - sub.Z.quantile(.25), sub.Z.quantile(.75) - sub.Z.median()])
ax.bar(range(len(labels)), vals, yerr=np.array(errs).T, capsize=3,
       color=[PAL.get(l, "#888") for l in labels], width=0.6)
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels([l.replace("GSE", "GSE\n") for l in labels], fontsize=6.6)
ax.set_ylabel("Gene-set Z score (median, IQR)")
ax.set_title("B  Enrichment centred at zero", loc="left", fontweight="bold", fontsize=8.2)

# --- C: 样本水平 PCA（新增，P0-3） ---
ax = fig.add_subplot(gs[0, 2])
X = pd.read_pickle(os.path.join(M, "GSE256068_expr.pkl"))
D = pd.read_csv(os.path.join(M, "GSE256068_design.csv"))
D = D[D["sample"].isin(X.columns)].reset_index(drop=True)
X = X[D["sample"].tolist()]
pc, var = pca(X)
for tis in ["Temporal", "Frontal", "Parietal", "Hippocampus"]:
    for g in ["Case", "Control"]:
        m = ((D.tissue == tis) & (D.group == g)).values
        if not m.any():
            continue
        ax.scatter(pc[m, 0], pc[m, 1], s=22, marker=TIS.get(tis, "o"),
                   facecolor=GRP[g], edgecolor="k", linewidth=0.35, alpha=0.85,
                   label="%s / %s (n=%d)" % (g, tis, int(m.sum())))
ax.set_xlabel("PC1 (%.1f%%)" % (var[0] * 100), fontsize=7.5)
ax.set_ylabel("PC2 (%.1f%%)" % (var[1] * 100), fontsize=7.5)
ax.axhline(0, color="grey", lw=0.5, ls=":"); ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.tick_params(labelsize=6.5)
ax.legend(fontsize=5.2, frameon=False, loc="upper left", handletextpad=0.25,
          labelspacing=0.25, borderpad=0.1)
ax.set_title("C  GSE256068 patient-level PCA:\ngroup is confounded with brain region",
             loc="left", fontweight="bold", fontsize=8.2)

fig.savefig(os.path.join(DST, "Fig2_deg_and_Z.png"), dpi=300, facecolor="white")
plt.close(fig)
from PIL import Image
im = Image.open(os.path.join(DST, "Fig2_deg_and_Z.png"))
im.convert("RGB").save(os.path.join(DST, "Fig2_deg_and_Z.tif"), "TIFF",
                       compression="tiff_lzw", dpi=(300, 300))
print("[Fig2] 三面板已重绘 %dx%d" % (im.width, im.height))

# ================= Supplementary Figure S1（四面板，重命名） =================
old = os.path.join(DST, "FigS_sample_structure.png")
new = os.path.join(DST, "FigS1_sample_structure.png")
if os.path.exists(old):
    os.replace(old, new)
    i2 = Image.open(new)
    i2.convert("RGB").save(os.path.join(DST, "FigS1_sample_structure.tif"), "TIFF",
                           compression="tiff_lzw", dpi=(300, 300))
    print("[FigS1] 重命名为 FigS1_sample_structure（%dx%d）" % (i2.width, i2.height))
for f in os.listdir(DST):
    print("   ", f)
