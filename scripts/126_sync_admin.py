#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 126: 同步全部 admin 文档（字数、核验日志、核查报告、清单、README、题名页）"""
import io, os, re

PKG = "."
ADM = os.path.join(PKG, "08_Admin")
fail = []


def edit(path, pairs, label):
    t = io.open(path, encoding="utf-8").read()
    for old, new in pairs:
        n = t.count(old)
        if n != 1:
            fail.append("[%s] %d 次: %s" % (label, n, old[:60])); continue
        t = t.replace(old, new)
    io.open(path, "w", encoding="utf-8", newline="").write(t)


# ---------- 1) Word_Count.md ----------
WC = os.path.join(ADM, "Word_Count.md")
t = io.open(WC, encoding="utf-8").read()
t = re.sub(r"\| Introduction \| [\d,]+ \|", "| Introduction | 541 |", t)
t = re.sub(r"\| Materials and methods \| [\d,]+ \|", "| Materials and methods | 1,400 |", t)
t = re.sub(r"\| Results \| [\d,]+ \|", "| Results | 1,832 |", t)
t = re.sub(r"\| Discussion \| [\d,]+ \|", "| Discussion | 1,058 |", t)
t = re.sub(r"\| Conclusions \| [\d,]+ \|", "| Conclusions | 154 |", t)
t = re.sub(r"\*\*4,833\*\*", "**4,914**", t)
t = re.sub(r"\*\*4,904\*\*", "**4,985**", t)
t = re.sub(r"\| 4,833 \|", "| 4,914 |", t)
t = re.sub(r"\| 4,904 \|", "| 4,985 |", t)
t = t.replace("4,833 words excl.", "4,914 words excl.").replace("4,904 incl.", "4,985 incl.")
io.open(WC, "w", encoding="utf-8", newline="").write(t)
print("Word_Count.md 已更新")

# ---------- 2) Verification_Log.md：加第 4 轮 ----------
VL = os.path.join(ADM, "Verification_Log.md")
v = io.open(VL, encoding="utf-8").read()
ADD = """

### Round 4 (figure and data audit)

A full audit of every figure, table and supporting file against the analysis outputs was performed (scripts 105-125). Two defects and three inconsistencies were found. Full detail in `Figure_Data_Audit.md`.

| # | Defect | Fix |
|---|---|---|
| 24 | **The tissue-matched result was misattributed.** The Results stated that a "tissue-matched analysis ... with tissue as a covariate" reduced the DEG rate "from 44.3% to 0.5% (-99%)". Re-running the models showed that 44.3% is itself the result of adding the tissue covariate, and that ~0.5% requires the six cell-type scores as well as tissue | Replaced with the verified model ladder: 59.6% -> 42.9% (tissue) -> 21.8% (cell scores) -> 0.4% (both) for the primary comparison, and 61.3% -> 44.3% -> 9.5% -> 1.0% restricted to non-epileptic cortex. The Discussion was corrected in the same way, and the new Table S12 records the full ladder |
| 25 | **Group and brain region are confounded in the largest cohort, and this was disclosed nowhere.** All 64 cases are temporal-lobe specimens whereas only 4 of 27 controls are region-matched (controls: 13 hippocampus, 8 frontal, 2 parietal, 4 temporal). Principal components 1 and 2 (36.0% of variance) separate brain regions, not disease groups | Disclosed in four places: Methods (a dedicated paragraph), Results, Discussion (as the first of four mechanisms), and Limitations (as the second, flagged "most importantly"). Added a "controls region-matched" column to Table 1 and Table S1, and added Figure 2C and Figure S1 showing the patient-level PCA |
| 26 | DEG counts were paired with the wrong denominators: the text read "11,923 DEGs (59.7% of 19,958 tested genes)", but 11,923 belongs to the 19,999-gene frame while 19,958 is the gene-set analysis frame (which gives 11,905) | Headline counts changed to 11,905 of 19,958 in the Abstract, Results and Key Points; the full two-frame statement moved to the Supplementary Methods and the inline note now states that each count carries its own denominator |
| 27 | The label-permutation instability for GSE213488 was not reported: the null 95th percentile differs six-fold between runs (3.25 with 100 permutations, 20.35 with 200), and Table S6 listed only DEG counts and P values | Table S6 extended with the null 95th percentile for both runs |
| 28 | No figure showed patient-level data, so the confounding in item 25 was invisible | Figure 2 gained panel C (patient-level PCA, group x region); Figure S1 added with four panels covering both cohorts and the adjustment ladder |

**Word-count consequence.** The disclosures and new Results text added ~430 words, taking the main text to 5,335. Four trim passes (scripts 122-125) removed ~350 words of duplicated wording without deleting any number, method detail or conclusion, moving two derivations (the analysis-set statement and the parametric-P rationale) into the Supplementary Methods where they already belonged. Final: **4,914 words excluding tables / 4,985 including**.

**Nothing was found to be fabricated.** Every value in Table 1, Tables S1-S11b and Figures 1-9 reproduced from the analysis outputs, with the single exception of the misattributed model in item 24.
"""
if "Round 4 (figure and data audit)" not in v:
    io.open(VL, "w", encoding="utf-8", newline="").write(v.rstrip() + "\n" + ADD)
    print("Verification_Log.md 已加第 4 轮")

# ---------- 3) Figure_Data_Audit.md：标注已修复 ----------
FA = os.path.join(ADM, "Figure_Data_Audit.md")
f = io.open(FA, encoding="utf-8").read()
f = f.replace("# 图表数据核查报告",
              "# 图表数据核查报告（已完成修改）\n\n"
              "> **状态：报告中的 P0/P1/P2 建议已全部执行完毕**（脚本 114–125）。\n"
              "> 执行后：正文 4,914 词（不含表）/ 4,985 词（含表）；新增 Figure 2C、Figure S1、Table S12；\n"
              "> Table 1 与 Table S1 增列脑区匹配度；TableS6 增列 null 95th percentile。\n"
              "> 修改记录见 `Verification_Log.md` 第 4 轮。")
io.open(FA, "w", encoding="utf-8", newline="").write(f)
print("Figure_Data_Audit.md 已标注状态")

# ---------- 4) Submission_Checklist.md ----------
edit(os.path.join(ADM, "Submission_Checklist.md"), [
    ("| 3 | Manuscripts exceeding the limit by \u22651000 words returned immediately | "
     "\u2705 4,833 words excl. tables / 4,904 incl. |",
     "| 3 | Manuscripts exceeding the limit by \u22651000 words returned immediately | "
     "\u2705 4,914 words excl. tables / 4,985 incl. |"),
    ("| All supplementary items cited in text and present as files | \u2705 Tables S1\u2013S11b "
     "+ 3 text files |",
     "| All supplementary items cited in text and present as files | \u2705 Tables S1\u2013S12, "
     "Figures S1, + 3 supplementary text files |"),
    ("| All figures cited in text in consecutive order | \u2705 Figures 1\u20139 |",
     "| All figures cited in text in consecutive order | \u2705 Figures 1\u20139 (Fig 2 now has "
     "panels A\u2013C) + Figure S1 |"),
    ("| Re-confirm requirements on the journal website before submitting, as policies can change. |",
     "| Re-confirm requirements on the journal website before submitting, as policies can change. |"),
], "checklist")
print("Submission_Checklist.md 已更新")

print()
if fail:
    print("!! 未命中:")
    for x in fail:
        print("   " + x)
else:
    print("全部 admin 文档已同步")
