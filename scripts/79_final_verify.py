#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 79: P1 完成后的全包最终核验"""
import os, re
import pandas as pd

PKG = "."
p = os.path.join(PKG, "03_Manuscript", "Manuscript_BIB.md")
t = open(p, encoding="utf-8").read()

print("=" * 88)
print(" P1 完成后的最终核验")
print("=" * 88)

# 1 引用一致性
body = t[:t.find("## References")]
refsec = re.search(r"## References(.*?)(?=## Figure legends)", t, re.S).group(1)
nums = sorted(set(int(m.group(1)) for m in re.finditer(r"^(\d+)\.\s", refsec, re.M)))
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
print("\n[1] 参考文献")
print("    条数 %d, 未被引用 %s, 悬空 %s" % (
    len(nums), [i for i in nums if i not in cited] or "无",
    [i for i in sorted(cited) if i not in nums] or "无"))

# 2 图表一致性
leg = t[t.find("## Figure legends"):t.find("## Supplementary")]
fc = set(int(m.group(1)) for m in re.finditer(r"Figure (\d+)", body))
fd = set(int(m.group(1)) for m in re.finditer(r"\*\*Figure (\d+)\.", leg))
print("\n[2] 图")
print("    正文引用 %s" % sorted(fc))
print("    图注定义 %s" % sorted(fd))
print("    未定义 %s / 未引用 %s" % (sorted(fc - fd) or "无", sorted(fd - fc) or "无"))

# 3 补充材料引用
sc = set(re.findall(r"Table S(\d+)", body))
sd = set(re.findall(r"Table S(\d+)", t[t.find("## Supplementary"):]))
print("\n[3] 补充表")
print("    正文引用 S%s" % sorted(sc, key=int))
print("    清单列出 S%s" % sorted(sd, key=int))
print("    未定义 %s" % (sorted(sc - sd, key=int) or "无"))

# 4 文件
print("\n[4] 补充表文件")
missing = [s for s in sorted(sc, key=int)
           if not any(f.startswith("TableS" + s + "_") for f in os.listdir(os.path.join(PKG, "06_Supplementary")))]
print("    缺失文件: %s" % (missing or "无"))

# 5 字数
def wc(s):
    return len(re.findall(r"\b[\w\u2019\-]+\b", s))


def sec(x, y):
    m = re.search(x + r"(.*?)(?=" + y + r"|\Z)", t, re.S | re.I)
    return m.group(1) if m else ""


parts = [sec(r"## 1\. Introduction", r"## 2\."), sec(r"## 2\. Materials and Methods", r"## 3\."),
         sec(r"## 3\. Results", r"## 4\."), sec(r"## 4\. Discussion", r"## Declarations")]
tot = sum(wc(x) for x in parts)
tb = sum(wc(" ".join(l for l in x.split("\n") if l.strip().startswith("|"))) for x in parts)
print("\n[5] 字数")
print("    正文(不含表格) %d / 5000 -> %s" % (tot - tb, "OK" if tot - tb <= 5000 else "超"))
print("    正文(含表格)   %d / 5000 -> %s" % (tot, "OK" if tot <= 5000 else "超 %d" % (tot - 5000)))

# 6 关键数字抽查
print("\n[6] 关键数字核对")
CHK = [
    ("1.6%", "1.6%"), ("576 篇", "576"), ("59.7%", "59.7%"), ("11,923", "11,923"),
    ("21/248", "21/248"), ("0/248", "0/248"), ("FDR < 10⁻⁴", "FDR < 10⁻⁴"),
    ("+31%", "+31%"), ("36%", "36%"), ("0 到 148", "0 to 148"),
    ("six-item", "six-item"), ("FDR = 0.24" if "0.24" in t else "0.24", "0.24"),
]
for lab, pat in CHK:
    print("    %-16s %s" % (lab, "✓" if pat in t else "✗ 未找到"))

# 7 包清单
print("\n[7] 投稿包")
n = 0; total = 0
for root, dirs, files in os.walk(PKG):
    for f in files:
        n += 1; total += os.path.getsize(os.path.join(root, f))
print("    文件数 %d, 合计 %.2f MB" % (n, total / 1e6))
for d in sorted(os.listdir(PKG)):
    fp = os.path.join(PKG, d)
    if os.path.isdir(fp):
        print("       %-22s %2d 个文件" % (d, len(os.listdir(fp))))
print("       %-22s %2d 个文件" % ("(根目录)", len([f for f in os.listdir(PKG) if os.path.isfile(os.path.join(PKG, f))])))
