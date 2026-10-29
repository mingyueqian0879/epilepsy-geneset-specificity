#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 59: P0 修改最终核验"""
import re, os

PKG = "."
p = os.path.join(PKG, "03_Manuscript", "Manuscript_BIB.md")
t = open(p, encoding="utf-8").read()

body = t[:t.find("## References")]
refsec = re.search(r"## References(.*?)(?=## Figure legends)", t, re.S).group(1)
nums = sorted(set(int(m.group(1)) for m in re.finditer(r"^(\d+)\.\s", refsec, re.M)))
cited = set()
for m in re.finditer(r"\[([\d,\u2013\-\s]+)\]", body):
    for part in m.group(1).split(","):
        part = part.strip().replace("\u2013", "-")
        if "-" in part:
            try:
                a, b = part.split("-")
                cited.update(range(int(a), int(b) + 1))
            except Exception:
                pass
        elif part.isdigit():
            cited.add(int(part))

print("=" * 78)
print(" P0 三项修改最终核验")
print("=" * 78)
print("  参考文献            %d 条" % len(nums))
print("  未被引用            %s" % ([i for i in nums if i not in cited] or "无"))
print("  悬空引用            %s" % ([i for i in sorted(cited) if i not in nums] or "无"))
print()
print("  P0-1  措辞改写")
print("        'do not exceed chance' 残留      %d 处 (应为 0)" % t.count("do not exceed chance"))
print("        'not evidence of its involvement'  %d 处 (保留 1 处, 属逻辑主张)" %
      t.count("not evidence of its involvement"))
print("        'could detect' / 'detectable'     %d 处" %
      len(re.findall(r"could detect|detectable|not estimable", t)))
print()
print("  P0-2  校准探针")
print("        'calibration probe(s)' 提及       %d 处" %
      (t.count("calibration probe") + t.count("calibration probes")))
print("        'rank' 报告方式                   %d 处" % len(re.findall(r"\brank(?:ed|s)?\b", t)))
print("        仍称 'positive control(s)' 于摘要/结果  %d 处" %
      len(re.findall(r"positive control", t)))
print()
print("  P0-3  MDE / 效力")
print("        'minimum detectable' 提及         %d 处" % len(re.findall(r"minimum detectable", t)))
print("        'MDE' 简称                        %d 处" % t.count("MDE"))
print("        'power' 提及                      %d 处" % len(re.findall(r"\bpower\b", t)))
print("        Table S7/S8/S9 交叉引用           %d 处" %
      len(re.findall(r"Table S[789]", t)))
print()
print("  报告标准条数        'six-item' %d 处 / 'five-item' %d 处" %
      (t.count("six-item"), t.count("five-item")))

# 词数
def wc(s):
    return len(re.findall(r"\b[\w\u2019\-]+\b", s))


def sec(a, b):
    m = re.search(a + r"(.*?)(?=" + b + r"|\Z)", t, re.S | re.I)
    return m.group(1) if m else ""


parts = [sec(r"## 1\. Introduction", r"## 2\."),
         sec(r"## 2\. Materials and Methods", r"## 3\."),
         sec(r"## 3\. Results", r"## 4\."),
         sec(r"## 4\. Discussion", r"## Declarations")]
tot = sum(wc(x) for x in parts)
tbl = sum(wc(" ".join(l for l in x.split("\n") if l.strip().startswith("|"))) for x in parts)
print()
print("  正文词数(含表格)    %d  -> %s" % (tot, "OK" if tot <= 5000 else "超 %d" % (tot - 5000)))
print("  正文词数(不含表格)  %d  -> %s" % (tot - tbl, "OK" if tot - tbl <= 5000 else "超"))

print()
print("=" * 78)
print(" 投稿包清点")
print("=" * 78)
n = 0
total = 0
for root, dirs, files in os.walk(PKG):
    if not files:
        continue
    rel = os.path.relpath(root, PKG)
    print("\n  %s/" % ("." if rel == "." else rel))
    for f in sorted(files):
        fp = os.path.join(root, f)
        sz = os.path.getsize(fp)
        n += 1
        total += sz
        print("     %-50s %8.1f KB" % (f, sz / 1024))
print("\n  文件数 = %d   合计 = %.2f MB" % (n, total / 1e6))
