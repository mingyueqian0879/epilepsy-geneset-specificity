#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 127: 同步 Title_Page / Cover_Letter / README，并重建封面信与题名页 DOCX"""
import io, os, subprocess, shutil

PKG = "."
fail = []


def edit(path, pairs, label):
    t = io.open(path, encoding="utf-8").read()
    for old, new in pairs:
        n = t.count(old)
        if n != 1:
            fail.append("[%s] %d 次: %s" % (label, n, old[:60])); continue
        t = t.replace(old, new)
    io.open(path, "w", encoding="utf-8", newline="").write(t)


# ---------- Title_Page ----------
edit(os.path.join(PKG, "02_Title_Page", "Title_Page.md"), [
    ("| Main text (Introduction\u2013Conclusions) | 4,833 words excl. tables (4,904 incl.) |",
     "| Main text (Introduction\u2013Conclusions) | 4,914 words excl. tables (4,985 incl.) |"),
    ("| Figures | 9 | \u2014 |", "| Figures | 9 + 1 supplementary (Figure S1) | \u2014 |"),
    ("| Supplementary items | 12 tables + 3 text files | \u2014 |",
     "| Supplementary items | 13 tables (S1\u2013S12) + Figure S1 + 3 text files | \u2014 |"),
], "title page")

# ---------- Cover Letter：统一 DEG 口径 + 补混杂披露 ----------
edit(os.path.join(PKG, "01_Cover_Letter", "Cover_Letter.md"), [
    ("we found that in the largest cohort 11,923 genes (59.7%) were differentially expressed, "
     "yet **none of 248 gene sets was enriched after FDR correction**.",
     "we found that in the largest cohort 11,905 genes (59.7%) were differentially expressed, "
     "yet **none of 248 gene sets was enriched after FDR correction**."),
    ("Label permutation confirmed the underlying differential-expression signal to be "
     "statistically real but **global** (empirical P = 0.005 in the largest cohort), and adjusting "
     "for cell composition removed 63% of DEGs.",
     "Label permutation confirmed the underlying differential-expression signal to be "
     "statistically real but **global** (empirical P = 0.005 in the largest cohort). We also "
     "report an important design limitation of that cohort openly: all 64 cases are temporal-lobe "
     "specimens whereas only 4 of 27 controls are region-matched, so part of the DEG rate reflects "
     "a regional contrast. The manuscript discloses this in the Methods, Results, Discussion and "
     "Limitations, adds region matching to Table 1, includes a patient-level principal-component "
     "panel as Figure 2C, and reports the adjustment ladder (59.6% \u2192 42.9% with tissue "
     "\u2192 21.8% with cell scores \u2192 0.4% with both) in Table S12."),
], "cover letter")

# ---------- README ----------
edit(os.path.join(PKG, "README.md"), [
    ("| Main text | 4,833 words excl. tables / 4,904 incl. (limit 2,000\u20135,000) |",
     "| Main text | 4,914 words excl. tables / 4,985 incl. (limit 2,000\u20135,000) |"),
    ("| Supplementary tables | 12 (S1\u2013S11b) + 3 text files |",
     "| Supplementary tables | 13 (S1\u2013S12) + Figure S1 + 3 text files |"),
    ("| Figures | 9, all 300 dpi |", "| Figures | 9 main (Fig 2 now has panels A\u2013C) + "
     "Figure S1, all 300 dpi |"),
], "readme")

if fail:
    print("!! 未命中:")
    for x in fail:
        print("   " + x)

# ---------- 重建封面信与题名页 DOCX（A4 双倍行距）----------
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

for rel in [r"01_Cover_Letter\Cover_Letter", r"02_Title_Page\Title_Page"]:
    md = os.path.join(PKG, rel + ".md")
    dx = os.path.join(PKG, rel + ".docx")
    r = subprocess.run(["pandoc", md, "-o", dx, "--from", "markdown", "--to", "docx"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("pandoc 失败:", rel); continue
    d = Document(dx)
    for s in d.sections:
        s.page_width, s.page_height = Cm(21.0), Cm(29.7)
        s.top_margin = s.bottom_margin = s.left_margin = s.right_margin = Cm(2.54)
    st = d.styles["Normal"]
    st.font.name = "Times New Roman"; st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    st.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    for p in d.paragraphs:
        p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        for run in p.runs:
            run.font.name = "Times New Roman"; run.font.size = Pt(12)
            rPr = run._element.get_or_add_rPr()
            rf = rPr.find(qn("w:rFonts"))
            if rf is None:
                rf = OxmlElement("w:rFonts"); rPr.insert(0, rf)
            for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
                rf.set(qn(a), "Times New Roman")
    d.save(dx)
    print("  重建:", os.path.basename(dx))

# ---------- 检查残留旧数字 ----------
print("\n=== 全包旧数字扫描 ===")
STALE = ["11,923 genes", "4,833", "4,904", "5036", "5038", "44.3% to 0.5%",
         "in a tissue-matched comparison 99%"]
for root, dirs, files in os.walk(PKG):
    for f in files:
        if not f.endswith((".md", ".csv")):
            continue
        p = os.path.join(root, f)
        s = io.open(p, encoding="utf-8", errors="ignore").read()
        for pat in STALE:
            if pat in s:
                print("  %s : %s" % (os.path.relpath(p, PKG), pat))
print("扫描完成")
