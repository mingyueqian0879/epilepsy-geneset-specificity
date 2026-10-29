#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 77: 生成 Table S10 (spiked-in) 与 S11 (方法对比), 并重建投稿包"""
import os, re
import pandas as pd

BASE = "."
OUT = os.path.join(BASE, "out")
PKG = "."
SUP = os.path.join(PKG, "06_Supplementary")
TAB = os.path.join(PKG, "05_Tables")

# ---------- Table S10: spiked-in ----------
sp = pd.read_csv(os.path.join(OUT, "spikein_simulation_corrected.csv"))
sp.columns = ["Spike magnitude (delta)", "Replicates", "Detection rate (%)",
              "Median Z", "Median observed DEGs", "Median expected DEGs",
              "Median achieved excess (%)", "False-positive rate at alpha=0.05 (%)"]
sp.to_csv(os.path.join(SUP, "TableS10_spikein_simulation.csv"), index=False)
print("[TableS10] 行数 =", len(sp))

# ---------- Table S11: 方法对比 ----------
mc = pd.read_csv(os.path.join(OUT, "method_comparison.csv"))
keep = ["dataset", "gene_set", "n", "Z_matched", "P_matched", "FDR_matched",
        "P_camera", "FDR_camera", "P_fry", "FDR_fry", "P_hyper", "FDR_hyper"]
mc = mc[[c for c in keep if c in mc.columns]].copy()
mc.columns = ["Cohort", "Gene set", "Genes detected", "Z (expression-matched)",
              "P (expression-matched)", "FDR (expression-matched)",
              "P (camera)", "FDR (camera)", "P (fry)", "FDR (fry)",
              "P (hypergeometric)", "FDR (hypergeometric)"][:len(mc.columns)]
mc = mc.sort_values(["Cohort", "FDR (expression-matched)"])
mc.to_csv(os.path.join(SUP, "TableS11_method_comparison.csv"), index=False)
print("[TableS11] 行数 =", len(mc))

# ---------- 汇总各方法显著数 ----------
rows = []
for ds in mc["Cohort"].unique():
    s = mc[mc["Cohort"] == ds]
    rows.append({"Cohort": ds, "Gene sets tested": len(s),
                 "Expression-matched (ours)": int((s["FDR (expression-matched)"] < 0.05).sum()),
                 "camera (competitive)": int((s["FDR (camera)"] < 0.05).sum()),
                 "fry (self-contained)": int((s["FDR (fry)"] < 0.05).sum()),
                 "Naive hypergeometric": int((s["FDR (hypergeometric)"] < 0.05).sum())})
sm = pd.DataFrame(rows)
sm.to_csv(os.path.join(SUP, "TableS11b_method_summary.csv"), index=False)
print("\n[TableS11b]")
print(sm.to_string(index=False))

# ---------- 更新稿件头部与题名页 ----------
p = os.path.join(PKG, "03_Manuscript", "Manuscript_BIB.md")
t = open(p, encoding="utf-8").read()


def wc(s):
    return len(re.findall(r"\b[\w\u2019\-]+\b", s))


def sec(x, y):
    m = re.search(x + r"(.*?)(?=" + y + r"|\Z)", t, re.S | re.I)
    return m.group(1) if m else ""


parts = [sec(r"## 1\. Introduction", r"## 2\."), sec(r"## 2\. Materials and Methods", r"## 3\."),
         sec(r"## 3\. Results", r"## 4\."), sec(r"## 4\. Discussion", r"## Declarations")]
tot = sum(wc(x) for x in parts)
tb = sum(wc(" ".join(l for l in x.split("\n") if l.strip().startswith("|"))) for x in parts)
noTbl = tot - tb

t = re.sub(r"\*\*Word count\*\*: main text [\d,]+ words[^;]*; abstract \d+ words; \d+ figures; 3 main-text tables; \d+ supplementary tables; 61 references",
           "**Word count**: main text %d words excluding tables (%d including inline tables), within the 2,000–5,000 limit for Problem solving protocols; abstract 296 words; 9 figures; 3 main-text tables; 12 supplementary items; 61 references" % (noTbl, tot), t)
open(p, "w", encoding="utf-8").write(t)
print("\n[稿件头部已更新] 正文 %d 词(不含表) / %d 词(含表)" % (noTbl, tot))

p2 = os.path.join(PKG, "02_Title_Page", "Title_Page.md")
t2 = open(p2, encoding="utf-8").read()
t2 = re.sub(r"\| Main text \(Introduction–Conclusions\) \| [^|]+\|",
            "| Main text (Introduction–Conclusions) | %d words excl. tables (%d incl.) |" % (noTbl, tot), t2)
t2 = re.sub(r"\| Abstract \| [^|]+\|", "| Abstract | 296 words |", t2)
t2 = re.sub(r"\| Figures \| \d+ \|", "| Figures | 9 |", t2)
t2 = re.sub(r"\| Supplementary tables \| \d+ \|", "| Supplementary items | 12 |", t2)
open(p2, "w", encoding="utf-8").write(t2)
print("[题名页已更新]")
