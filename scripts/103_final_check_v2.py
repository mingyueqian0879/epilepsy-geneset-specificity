#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 103: 新投稿包终检"""
import io, os, re
from docx import Document
from docx.enum.text import WD_LINE_SPACING

NEW = "."
MAN = os.path.join(NEW, "03_Manuscript", "Manuscript_BIBformat.docx")
MD = os.path.join(NEW, "03_Manuscript", "Manuscript_BIBformat.md")
fail = []

print("=" * 92)
print(" BIB 体例投稿包终检")
print("=" * 92)

# 1 文件清单
n_file = n_mb = 0
for root, dirs, files in os.walk(NEW):
    for f in files:
        n_file += 1; n_mb += os.path.getsize(os.path.join(root, f))
print("\n[1] 文件清单")
print("    共 %d 个文件, %.2f MB" % (n_file, n_mb / 1e6))
for d in sorted(os.listdir(NEW)):
    p = os.path.join(NEW, d)
    if os.path.isdir(p):
        print("      %-20s %2d" % (d, len(os.listdir(p))))
    else:
        print("      %-20s  (根文件)" % d)

# 2 Word 稿件格式
print("\n[2] Word 稿件格式（期刊要求：Word、双倍行距、A4）")
d = Document(MAN)
s0 = d.sections[0]
a4 = abs(s0.page_width.cm - 21.0) < 0.1 and abs(s0.page_height.cm - 29.7) < 0.1
dbl = all(p.paragraph_format.line_spacing_rule == WD_LINE_SPACING.DOUBLE for p in d.paragraphs)
print("    A4 页面: %s (%.1f x %.1f cm)" % ("是" if a4 else "否", s0.page_width.cm, s0.page_height.cm))
print("    双倍行距: %s" % ("是" if dbl else "否"))
print("    字体: %s %.0fpt" % (d.styles["Normal"].font.name, d.styles["Normal"].font.size.pt))
print("    页边距: %.2f cm" % s0.top_margin.cm)
if not (a4 and dbl):
    fail.append("Word 格式不符")

# 3 结构
t = io.open(MD, encoding="utf-8").read()
tops = re.findall(r"^## (.+)$", t, re.M)
expect = ["Abstract", "Introduction", "Materials and methods", "Results", "Discussion",
          "Conclusions", "Supplementary material", "Acknowledgements", "Author contributions",
          "Conflict of interest", "Funding", "Data availability",
          "Ethics approval and consent to participate", "Consent for publication",
          "Use of artificial intelligence", "References", "Figure legends"]
print("\n[3] 章节结构")
print("    一级章节 %d 个，与 BIB 体例一致: %s" % (len(tops), "是" if tops == expect else "否"))
if tops != expect:
    fail.append("章节顺序不符")
print("    带编号的标题:", re.findall(r"^#{2,3} \d+\.", t, re.M) or "无")
print("    残留 § 交叉引用:", sorted(set(re.findall(r"\u00a7[\d.]+", t))) or "无")
if re.findall(r"\u00a7[\d.]+", t):
    fail.append("残留编号引用")

# 4 参考文献
body = t[:t.find("## References")]
refsec = re.search(r"## References(.*?)(?=## Figure legends)", t, re.S).group(1)
nums = sorted(int(m.group(1)) for m in re.finditer(r"^(\d+)\\?\.\s", refsec, re.M))
cited = set()
for m in re.finditer(r"\[([\d,\u2013\-\s]+)\]", body):
    for part in m.group(1).split(","):
        part = part.strip().replace("\u2013", "-")
        if "-" in part:
            try:
                x, y = part.split("-"); cited.update(range(int(x), int(y) + 1))
            except Exception:
                pass
        elif part.isdigit():
            cited.add(int(part))
ndoi = len(re.findall(r"10\.\d{4,}/", refsec))
nissue = len([l for l in refsec.split("\n")
              if re.match(r"^\d+\\?\.", l.strip()) and re.search(r";\d+\(\d+\)", l)])
print("\n[4] 参考文献（Oxford SCIMED）")
print("    条数 %d，编号连续: %s" % (len(nums), "是" if nums == list(range(1, len(nums) + 1)) else "否"))
print("    未被引用: %s | 悬空: %s" % ([i for i in nums if i not in cited] or "无",
                                       [i for i in sorted(cited) if i not in nums] or "无"))
print("    带 DOI: %d/%d | 带期号(应为0): %d" % (ndoi, len(nums), nissue))
if [i for i in nums if i not in cited] or nissue:
    fail.append("参考文献问题")

# 5 图与补充表
fc = sorted(set(int(m.group(1)) for m in re.finditer(r"Figure (\d+)", body)))
fd = sorted(set(int(m.group(1)) for m in re.finditer(r"\*\*Figure (\d+)\.", t)))
sc = set(re.findall(r"Table S(\d+[ab]?)", body))
sdir = os.path.join(NEW, "06_Supplementary")
miss = [s for s in sc if not any(f.startswith("TableS" + s + "_") for f in os.listdir(sdir))]
print("\n[5] 图 / 补充材料")
print("    图 %s，图注 %s，未定义: %s" % (fc, fd, sorted(set(fc) - set(fd)) or "无"))
print("    正文引用补充表 %d 项，缺文件: %s" % (len(sc), miss or "无"))
if sorted(set(fc) - set(fd)) or miss:
    fail.append("图表问题")

# 6 字数
def wc(s):
    return len(re.findall(r"\b[\w\u2019\-]+\b", s))


def sec(x, y):
    m = re.search(x + r"(.*?)(?=" + y + r"|\Z)", t, re.S | re.I)
    return m.group(1) if m else ""


conc = sec(r"## Conclusions", r"## Supplementary material")
k = conc.find("> **Key Points**")
parts = [sec(r"## Introduction", r"## Materials and methods"),
         sec(r"## Materials and methods", r"## Results"),
         sec(r"## Results", r"## Discussion"),
         sec(r"## Discussion", r"## Conclusions"), conc[:k] if k > 0 else conc]
tot = sum(wc(x) for x in parts)
tb = sum(wc(" ".join(l for l in x.split("\n") if l.strip().startswith("|"))) for x in parts)
print("\n[6] 字数")
print("    正文不含表格 %d / 5000 -> %s" % (tot - tb, "OK" if tot - tb <= 5000 else "超"))
print("    正文含表格   %d / 5000 -> %s" % (tot, "OK" if tot <= 5000 else "超"))
print("    Key Points %d 词" % wc(conc[k:] if k > 0 else ""))
if tot > 5000:
    fail.append("字数超限")

# 7 关键数字
print("\n[7] 关键数字")
KEY = [("1.6%", "1.6%"), ("576", "576"), ("59.7%", "59.7%"), ("11,923", "11,923"),
       ("21/248", "21/248"), ("0/248", "0/248"), ("36%", "36%"), ("0 to 148", "0 to 148"),
       ("six-item", "six-item"), ("0.24", "0.24"), ("+31", "31%"), ("spiked", "piked-in")]
for lab, pat in KEY:
    ok = pat in t
    print("    %-14s %s" % (lab, "OK" if ok else "** 未找到 **"))
    if not ok:
        fail.append("缺关键数字 " + lab)

# 8 admin 文档无过期数字
print("\n[8] admin 文档过期数字扫描")
STALE = ["5036", "5038", "3,794", "18/18", "Fig1\u2026Fig8", "S1\u2013S5", "Figures 1\u20138",
         "278 word", "6/248", "r = \u22120.012", "\u00a72.", "\u00a73."]
hit = []
for root, dirs, files in os.walk(NEW):
    for f in files:
        if not f.endswith(".md"):
            continue
        p = os.path.join(root, f)
        s = io.open(p, encoding="utf-8", errors="ignore").read()
        for pat in STALE:
            if pat in s:
                hit.append("%s: %s" % (os.path.relpath(p, NEW), pat))
print("    命中:", hit or "无")

print("\n" + "=" * 92)
print(" 结论:", "全部通过" if not fail else "需处理 -> %s" % fail)
print("=" * 92)
