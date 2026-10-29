#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 127: 修复 Fig 2 与 Figure S1 的三处版面/数值缺陷，并自检。

修复项
------
(1) Fig 2A  x 轴第三行 "region matched: 4/27" 与相邻刻度标签水平重叠
    （实测每对相邻标签重叠 38-44 px）→ 改为 "region 4/27"，字号 6.4 → 6.0
(2) Fig 2C  图例（5 项）压在散点上（实测盖住 4 个样本点）
    → 2 列 + y 轴留白 18%，实测遮挡 0
(3) Fig S1A 同一问题（实测盖住 5 个样本点）→ 2 列 + y 轴留白 20%
(4) Fig S1D 阶梯数值错误：图中为 59.6/42.9/**17.0**/0.4，
    但 limma 复跑（脚本 115）与 Table S12 均为 59.62/42.88/**21.84**/0.40，
    正文亦写 21.8%。→ 改为直接从 TableS12_model_ladder.csv 读取，杜绝再次漂移。

自检：图例遮挡按「散点圆心 + 标记半径」判定（比只判圆心更严格）；
      x 轴标签按渲染宽度 vs 刻度间距判定。任一不过即抛错。
"""
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

BASE = "."
OUT = os.path.join(BASE, "out")
M = os.path.join(OUT, "multi")
DST = "figures"
S12 = "results/supplementary/TableS12_model_ladder.csv"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.linewidth": 0.8, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300,
                     "savefig.bbox": "tight", "savefig.pad_inches": 0.06})
PAL = {"GSE256068": "#2b6cb0", "GSE213488": "#c05621",
       "GSE128300": "#2f855a", "GSE272142": "#6b46c1"}
GRP = {"Case": "#c53030", "Control": "#2b6cb0"}
TIS = {"Temporal": "o", "Frontal": "s", "Parietal": "^", "Hippocampus": "D",
       "grey_matter": "o", "cortex": "o", "dentate_gyrus": "o"}


def pca(X):
    A = X.values.astype(float)
    A = A - A.mean(axis=1, keepdims=True)
    sd = A.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    A = np.nan_to_num(A / sd)
    _, S_, Vt = np.linalg.svd(A, full_matrices=False)
    return (Vt.T * S_)[:, :2], (S_ ** 2) / (S_ ** 2).sum()


def load_256068():
    X = pd.read_pickle(os.path.join(M, "GSE256068_expr.pkl"))
    D = pd.read_csv(os.path.join(M, "GSE256068_design.csv"))
    D = D[D["sample"].isin(X.columns)].reset_index(drop=True)
    return X[D["sample"].tolist()], D


def scatter_groups(ax, pc, D, order, size, fs_labels=True):
    """按 (tissue, group) 逐组打点并返回 (labels, masks)。"""
    labs, masks = [], []
    for tis in order:
        for g in ["Case", "Control"]:
            m = ((D.tissue == tis) & (D.group == g)).values
            if not m.any():
                continue
            labs.append("%s / %s (n=%d)" % (g, tis, int(m.sum())))
            masks.append(m)
            ax.scatter(pc[m, 0], pc[m, 1], s=size, marker=TIS.get(tis, "o"),
                       facecolor=GRP[g], edgecolor="k", linewidth=0.35,
                       alpha=0.85, label=labs[-1])
    return labs, masks


def check_legend(ax, leg, pc, size, tag):
    """严格遮挡判定：散点圆心 + 标记半径 落入图例矩形即算遮挡。"""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = leg.get_window_extent(r)
    radius_px = (np.sqrt(size) / 2.0) * (fig.dpi / 72.0)
    disp = ax.transData.transform(np.asarray(pc)[:, :2])
    inside = ((disp[:, 0] >= bb.x0 - radius_px) & (disp[:, 0] <= bb.x1 + radius_px)
              & (disp[:, 1] >= bb.y0 - radius_px) & (disp[:, 1] <= bb.y1 + radius_px))
    inv = ax.transAxes.inverted()
    (x0, y0), (x1, y1) = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    n = int(inside.sum())
    print("    [%s] occluded=%d (marker-radius margin %.1f px)  bbox x=[%.2f,%.2f]"
          % (tag, n, radius_px, x0, x1))
    if n:
        raise SystemExit("FAIL: %s legend still occludes %d point(s)" % (tag, n))
    if x1 > 0.995:
        raise SystemExit("FAIL: %s legend crosses right spine (x1=%.3f)" % (tag, x1))
    return n


def check_xticklabels(ax, tag):
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bbs = [t.get_window_extent(r) for t in ax.get_xticklabels()]
    worst = max(bbs[i].x1 - bbs[i + 1].x0 for i in range(len(bbs) - 1))
    print("    [%s] worst adjacent-label overlap = %+.1f px" % (tag, worst))
    if worst > 0:
        raise SystemExit("FAIL: %s x tick labels overlap by %.1f px" % (tag, worst))
    return worst


# =====================================================================
# 读取权威阶梯值（Table S12），不再硬编码
# =====================================================================
t12 = pd.read_csv(S12)
prim = t12[t12["Comparison"].str.startswith("GSE256068: 64 cases vs all 27 controls")]
prim = prim.sort_values("Model")
LADDER = list(prim["DEG rate (%)"])
assert len(LADDER) == 4, LADDER
print("Step 127 — 修复 Fig 2 / Figure S1")
print("  Table S12 阶梯（64 cases vs 27 controls）: " +
      " -> ".join("%.2f%%" % v for v in LADDER))

# =====================================================================
# Fig 2（三面板）
# =====================================================================
print("\n[Fig 2]")
fig = plt.figure(figsize=(11.0, 3.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.85, 1.15], wspace=0.34)

S = pd.read_csv(os.path.join(OUT, "dataset_summary.csv"))
R = pd.read_csv(os.path.join(OUT, "specificity_matrix_final.csv"))
MATCH = {"GSE256068": "4/27", "GSE213488": "4/4",
         "GSE128300": "8/8", "GSE272142": "6/6"}

# --- A ---
axA = fig.add_subplot(gs[0, 0])
ds = S.sort_values("deg_rate", ascending=False)
axA.bar(range(len(ds)), ds.deg_rate,
        color=[PAL.get(d, "#888") for d in ds.dataset], width=0.62)
axA.set_xticks(range(len(ds)))
axA.set_xticklabels(["%s\n(n=%d/%d)\nregion %s"
                     % (d, c, k, MATCH.get(d, "?"))
                     for d, c, k in zip(ds.dataset, ds.n_case, ds.n_ctrl)],
                    fontsize=6.0)
axA.set_ylabel("DEG rate at FDR<0.05 (%)")
axA.axhspan(5, 20, color="gray", alpha=0.18, zorder=0)
axA.text(len(ds) - 0.5, 12.5, "typical\n5-20%", fontsize=7, ha="right",
         color="#4a5568")
for i, v in enumerate(ds.deg_rate):
    axA.text(i, v + 1.2, "%.1f%%" % v, ha="center", fontsize=8)
axA.set_ylim(0, 70)
axA.set_title("A  Same clinical question, DEG rate spans 0-60%", loc="left",
              fontweight="bold", fontsize=8.2)

# --- B ---
axB = fig.add_subplot(gs[0, 1])
labels, vals, errs = [], [], []
for d in ds.dataset:
    sub = R[R.dataset == d].dropna(subset=["Z"])
    if not len(sub):
        continue
    labels.append(d)
    vals.append(sub.Z.median())
    errs.append([sub.Z.median() - sub.Z.quantile(.25),
                 sub.Z.quantile(.75) - sub.Z.median()])
axB.bar(range(len(labels)), vals, yerr=np.array(errs).T, capsize=3,
        color=[PAL.get(l, "#888") for l in labels], width=0.6)
axB.axhline(0, color="k", lw=0.8)
axB.set_xticks(range(len(labels)))
axB.set_xticklabels([l.replace("GSE", "GSE\n") for l in labels], fontsize=6.6)
axB.set_ylabel("Gene-set Z score (median, IQR)")
axB.set_title("B  Enrichment centred at zero", loc="left", fontweight="bold",
              fontsize=8.2)

# --- C ---
axC = fig.add_subplot(gs[0, 2])
X, D = load_256068()
pc, var = pca(X)
labsC, _ = scatter_groups(axC, pc, D,
                          ["Temporal", "Frontal", "Parietal", "Hippocampus"], 22)
axC.set_xlabel("PC1 (%.1f%%)" % (var[0] * 100), fontsize=7.5)
axC.set_ylabel("PC2 (%.1f%%)" % (var[1] * 100), fontsize=7.5)
axC.axhline(0, color="grey", lw=0.5, ls=":")
axC.axvline(0, color="grey", lw=0.5, ls=":")
axC.tick_params(labelsize=6.5)
lo, hi = axC.get_ylim()
axC.set_ylim(lo, hi + 0.18 * (hi - lo))          # 顶部留白给图例
legC = axC.legend(fontsize=5.8, ncol=2, frameon=False, loc="upper left",
                  handlelength=1.2, handletextpad=0.25, columnspacing=0.4,
                  labelspacing=0.3, borderpad=0.15)
axC.set_title("C  GSE256068 patient-level PCA:\ngroup is confounded with brain region",
              loc="left", fontweight="bold", fontsize=8.2)

check_xticklabels(axA, "Fig 2A tick labels")
check_legend(axC, legC, pc, 22, "Fig 2C legend")

fig.savefig(os.path.join(DST, "Fig2_deg_and_Z.png"), dpi=300, facecolor="white")
plt.close(fig)
im = Image.open(os.path.join(DST, "Fig2_deg_and_Z.png"))
im.convert("RGB").save(os.path.join(DST, "Fig2_deg_and_Z.tif"), "TIFF",
                       compression="tiff_lzw", dpi=(300, 300))
print("    Fig2 已重绘 %dx%d" % (im.width, im.height))

# =====================================================================
# Figure S1（四面板）
# =====================================================================
print("\n[Fig S1]")
fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0))

# --- A ---
ax = axes[0, 0]
labsA, _ = scatter_groups(ax, pc, D,
                          ["Temporal", "Frontal", "Parietal", "Hippocampus"], 26)
ax.set_xlabel("PC1 (%.1f%%)" % (var[0] * 100))
ax.set_ylabel("PC2 (%.1f%%)" % (var[1] * 100))
ax.axhline(0, color="grey", lw=0.5, ls=":")
ax.axvline(0, color="grey", lw=0.5, ls=":")
lo, hi = ax.get_ylim()
ax.set_ylim(lo, hi + 0.20 * (hi - lo))
legA = ax.legend(fontsize=6.3, ncol=2, frameon=False, loc="upper left",
                 handlelength=1.2, handletextpad=0.3, columnspacing=0.5,
                 labelspacing=0.3)
ax.set_title("A  GSE256068: group is confounded with brain region",
             loc="left", fontweight="bold", fontsize=8.8)

# --- B ---
ax = axes[0, 1]
ct = pd.crosstab(D.tissue, D.group)
ct = ct.reindex(columns=[c for c in ["Case", "Control"] if c in ct.columns])
ct.plot(kind="barh", stacked=True, ax=ax,
        color=[GRP[c] for c in ct.columns], width=0.62, edgecolor="k",
        linewidth=0.4)
for i, tis in enumerate(ct.index):
    tot = ct.loc[tis].sum()
    ax.text(tot + 1.2, i, "n=%d" % int(tot), va="center", fontsize=7)
ax.set_xlabel("Samples")
ax.set_ylabel("")
ax.legend(fontsize=7, frameon=False, loc="lower right")
ax.set_title("B  Only 4 of 27 controls match the cases' region",
             loc="left", fontweight="bold", fontsize=8.8)

# --- C ---
ax = axes[1, 0]
X2 = pd.read_pickle(os.path.join(M, "GSE213488_expr.pkl"))
D2 = pd.read_csv(os.path.join(M, "GSE213488_design.csv"))
D2 = D2[D2["sample"].isin(X2.columns)].reset_index(drop=True)
X2 = X2[D2["sample"].tolist()]
pc2, var2 = pca(X2)
for g in ["Case", "Control"]:
    m = (D2.group == g).values
    ax.scatter(pc2[m, 0], pc2[m, 1], s=34, marker="o", facecolor=GRP[g],
               edgecolor="k", linewidth=0.4, alpha=0.85,
               label="%s (n=%d)" % (g, m.sum()))
ax.set_xlabel("PC1 (%.1f%%)" % (var2[0] * 100))
ax.set_ylabel("PC2 (%.1f%%)" % (var2[1] * 100))
ax.axhline(0, color="grey", lw=0.5, ls=":")
ax.axvline(0, color="grey", lw=0.5, ls=":")
ax.legend(fontsize=7, frameon=False)
ax.set_title("C  GSE213488: single tissue, no such confounding",
             loc="left", fontweight="bold", fontsize=8.8)

# --- D ---
ax = axes[1, 1]
labels = ["all controls\n(no covariate)", "+ tissue\ncovariate",
          "+ 6 cell scores", "+ tissue\n+ 6 cell scores"]
cols = ["#a0aec0", "#90cdf4", "#63b3ed", "#2b6cb0"]
ax.bar(range(4), LADDER, color=cols, width=0.62, edgecolor="k", linewidth=0.4)
for i, v in enumerate(LADDER):
    ax.text(i, v + 1.5, "%.1f%%" % v, ha="center", fontsize=7.5)
ax.set_xticks(range(4))
ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("Genome-wide DEG rate (%)")
ax.set_ylim(0, 70)
ax.set_title("D  GSE256068: adjustment removes almost all signal",
             loc="left", fontweight="bold", fontsize=8.8)

fig.suptitle("Sample-level structure of the cohorts "
             "(patient-level evidence absent from Figures 2 and 3)",
             fontsize=10.2, y=1.02, x=0.02, ha="left")
fig.tight_layout()
fig.canvas.draw()
check_legend(axes[0, 0], legA, pc, 26, "Fig S1A legend")
print("    ladder printed in panel D: " +
      " / ".join("%.1f%%" % v for v in LADDER))

out_png = os.path.join(DST, "FigS1_sample_structure.png")
fig.savefig(out_png, dpi=300, facecolor="white")
plt.close(fig)
im2 = Image.open(out_png)
im2.convert("RGB").save(os.path.join(DST, "FigS1_sample_structure.tif"), "TIFF",
                        compression="tiff_lzw", dpi=(300, 300))
print("    FigS1 已重绘 %dx%d" % (im2.width, im2.height))

print("\n全部自检通过：标签无重叠、图例零遮挡、阶梯值取自 Table S12。")
