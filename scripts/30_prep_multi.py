#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 30: 统一预处理 4 个独立数据集 -> 标准化 (表达矩阵 + 分组表)
  GSE256068 : log2 TPM, TLE-HS vs Control
  GSE213488 : raw counts(grey matter), FCD vs Control(autopsy)
  GSE128300 : raw counts, Epilepsy(FCD+HS) vs Non-epilepsy
  GSE272142 : raw counts, TLE vs control
"""
import gzip, os, re, json
import numpy as np
import pandas as pd

BASE = "."
DL = os.path.join(BASE, "data", "multi")
OUT = os.path.join(BASE, "out", "multi")
os.makedirs(OUT, exist_ok=True)

meta = json.load(open(os.path.join(BASE, "out", "multi_design_raw.json"), encoding="utf-8"))


def save(name, expr, design):
    """expr: genes x samples (DataFrame), design: DataFrame(gsm/sample, group, dataset)"""
    e = expr.copy()
    e.index = [str(i).split(".")[0] for i in e.index]
    e = e[~e.index.duplicated(keep="first")]
    e.to_pickle(os.path.join(OUT, name + "_expr.pkl"))
    design.to_csv(os.path.join(OUT, name + "_design.csv"), index=False, encoding="utf-8-sig")
    print("   [saved] %s  expr=%s  design=%s" % (name, e.shape, design.shape))


# ---------------- 1. GSE256068 ----------------
print("\n[1] GSE256068")
expr = pd.read_pickle(os.path.join(BASE, "out", "expr_log2tpm.pkl"))
d = pd.read_csv(os.path.join(BASE, "out", "design_aligned.csv"), encoding="utf-8-sig")
assert "matrix_col" in d.columns, d.columns.tolist()
keep_sets = ["TLE-HS", "ControlCortex", "ControlHippocampus"]
sub = d[d["disease"].isin(keep_sets) & d["matrix_col"].isin(expr.columns)].copy()
sub = sub.drop_duplicates("matrix_col")
sub["group"] = np.where(sub["disease"] == "TLE-HS", "Case", "Control")
e1 = expr[sub["matrix_col"].tolist()]
save("GSE256068", e1,
     pd.DataFrame({"sample": sub["matrix_col"], "group": sub["group"],
                   "detail": sub["disease"], "tissue": sub["tissue"],
                   "dataset": "GSE256068"}))

# ---------------- 2. GSE213488 (grey matter) ----------------
print("\n[2] GSE213488 (grey matter)")
f = os.path.join(DL, "GSE213488_grey_matter_count.csv.gz")
g = pd.read_csv(f, index_col=0)
print("   counts shape =", g.shape)
# 用 Sample_description (s1_gm...) 匹配列名 s1, s3...
mg = meta["GSE213488"]
desc2dis = {}
for de, st in zip(mg["desc"], [None] * len(mg["desc"])):
    pass
# 从 series matrix 重新取 characteristics 与 description 的对应
tags = {}
raw = None
import urllib.request
u = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213488/matrix/"
     "GSE213488_series_matrix.txt.gz")
raw = gzip.decompress(urllib.request.urlopen(
    urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"}), timeout=180).read()
).decode("utf-8", "replace")
for line in raw.split("\n"):
    if line.startswith("!series_matrix_table_begin"):
        break
    if line.startswith("!Sample_"):
        p = line.rstrip().split("\t")
        tags.setdefault(p[0], []).append([x.strip().strip('"') for x in p[1:]])

desc = None
dis = None
for k in tags:
    if "description" in k:
        desc = tags[k][0]
    if "characteristics" in k:
        for v in tags[k]:
            if v and v[0].lower().startswith("disease state"):
                dis = [x.split(":", 1)[1].strip() for x in v]
m = pd.DataFrame({"desc": desc, "disease": dis}).dropna()
m = m[m["desc"].str.endswith("_gm")].copy()
m["key"] = m["desc"].str.replace("_gm", "", regex=False)
print("   grey-matter samples =", len(m))
print(m["disease"].value_counts().to_string())
cols = [c for c in g.columns if c in set(m["key"])]
print("   匹配列 =", len(cols), "/", g.shape[1])
gm = m[m["key"].isin(cols)].set_index("key").loc[cols]
gm["group"] = np.where(gm["disease"].str.contains("Control", case=False), "Control", "Case")
save("GSE213488", g[cols],
     pd.DataFrame({"sample": cols, "group": gm["group"].values,
                   "detail": gm["disease"].values, "tissue": "grey_matter",
                   "dataset": "GSE213488"}))

# ---------------- 3. GSE128300 ----------------
print("\n[3] GSE128300")
f = os.path.join(DL, "GSE128300_RNAstrandedRev_simp_2019.mx.txt.gz")
g3 = pd.read_csv(f, sep="\t", index_col=0)
print("   shape =", g3.shape)
mg3 = meta["GSE128300"]
# title = KK1-RNA ; 列名 = KK1-MRNA ; 去掉后缀对齐
t2d = {}
for t, dd in zip(mg3["title"], [None] * len(mg3["title"])):
    pass
u3 = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/GSE128300/matrix/"
      "GSE128300_series_matrix.txt.gz")
raw3 = gzip.decompress(urllib.request.urlopen(
    urllib.request.Request(u3, headers={"User-Agent": "Mozilla/5.0"}), timeout=180).read()
).decode("utf-8", "replace")
tg3 = {}
for line in raw3.split("\n"):
    if line.startswith("!series_matrix_table_begin"):
        break
    if line.startswith("!Sample_"):
        p = line.rstrip().split("\t")
        tg3.setdefault(p[0], []).append([x.strip().strip('"') for x in p[1:]])
ttl = None; dss = None
for k in tg3:
    if "title" in k:
        ttl = tg3[k][0]
    if "characteristics" in k:
        for v in tg3[k]:
            if v and v[0].lower().startswith("disease:"):
                dss = [x.split(":", 1)[1].strip() for x in v]
m3 = pd.DataFrame({"title": ttl, "disease": dss})
m3["key"] = m3["title"].str.replace("-RNA", "", regex=False).str.strip()
colmap = {}
for c in g3.columns:
    k = c.replace("-MRNA", "").strip()
    if k in set(m3["key"]):
        colmap[c] = k
print("   匹配列 =", len(colmap), "/", g3.shape[1])
sel = list(colmap.keys())
info = m3.set_index("key").loc[[colmap[c] for c in sel]]
info["group"] = np.where(info["disease"].str.contains("Non-epilepsy", case=False),
                         "Control", "Case")
save("GSE128300", g3[sel],
     pd.DataFrame({"sample": sel, "group": info["group"].values,
                   "detail": info["disease"].values, "tissue": "cortex",
                   "dataset": "GSE128300"}))

# ---------------- 4. GSE272142 ----------------
print("\n[4] GSE272142")
f = os.path.join(DL, "GSE272142_genes.rawcount.anno.xlsx")
g4 = pd.read_excel(f)
g4 = g4.set_index("#id")
ann_cols = ["Entrez", "Uniprot", "Position_gene", "Symbol", "Description"]
sample_cols = [c for c in g4.columns if c not in ann_cols]
counts4 = g4[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
print("   shape =", counts4.shape, " samples =", len(sample_cols))
m4 = meta["GSE272142"]  # characteristics 通道2 = disease state
u4 = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE272nnn/GSE272142/matrix/"
      "GSE272142_series_matrix.txt.gz")
raw4 = gzip.decompress(urllib.request.urlopen(
    urllib.request.Request(u4, headers={"User-Agent": "Mozilla/5.0"}), timeout=180).read()
).decode("utf-8", "replace")
tg4 = {}
for line in raw4.split("\n"):
    if line.startswith("!series_matrix_table_begin"):
        break
    if line.startswith("!Sample_"):
        p = line.rstrip().split("\t")
        tg4.setdefault(p[0], []).append([x.strip().strip('"') for x in p[1:]])
ttl4 = tg4.get("!Sample_title", [[""]])[0]
ds4 = None
for k in tg4:
    if "characteristics" in k:
        for v in tg4[k]:
            if v and v[0].lower().startswith("disease state"):
                ds4 = [x.split(":", 1)[1].strip() for x in v]
m4 = pd.DataFrame({"title": ttl4, "disease": ds4})
print(m4["disease"].value_counts().to_string())
print("   矩阵列 =", sample_cols[:8])
print("   元数据 title =", list(m4["title"])[:8])
# 列名与 title 直接匹配
sel4 = [c for c in sample_cols if c in set(m4["title"])]
print("   直接匹配 =", len(sel4), "/", len(sample_cols))
if len(sel4) < 6:
    # 按顺序匹配
    sel4 = sample_cols
    m4i = m4.copy()
    print("   -> 按顺序匹配")
    info4 = pd.DataFrame({"sample": sel4, "disease": m4["disease"].values[:len(sel4)]})
else:
    info4 = m4.set_index("title").loc[sel4].reset_index()
    info4.columns = ["sample", "disease"]
info4["group"] = np.where(info4["disease"].str.contains("epilepsy", case=False),
                          "Case", "Control")
save("GSE272142", counts4[info4["sample"].tolist()],
     pd.DataFrame({"sample": info4["sample"], "group": info4["group"],
                   "detail": info4["disease"], "tissue": "dentate_gyrus",
                   "dataset": "GSE272142"}))

print("\n" + "=" * 78)
print(" 汇总")
print("=" * 78)
for n in ["GSE256068", "GSE213488", "GSE128300", "GSE272142"]:
    dp = os.path.join(OUT, n + "_design.csv")
    if os.path.exists(dp):
        d = pd.read_csv(dp, encoding="utf-8-sig")
        print("  %-12s Case=%-4d Control=%-4d  (%s)" % (
            n, (d["group"] == "Case").sum(), (d["group"] == "Control").sum(),
            d["tissue"].iloc[0]))
