#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 120: 把 TableS12 接入正文 + 重建 Word 投稿文件"""
import io, os, re, subprocess, sys

PKG = "."
MS = os.path.join(PKG, "03_Manuscript", "Manuscript_BIBformat.md")
t = io.open(MS, encoding="utf-8").read()
fail = []


def rep(old, new, tag):
    global t
    n = t.count(old)
    if n != 1:
        fail.append("[%s] %d 次: %s" % (tag, n, old[:70])); return
    t = t.replace(old, new)


# 补充材料清单加 Table S12
rep("| **Table S11b** | Method comparison summary: number of significant sets per test and cohort |",
    "| **Table S11b** | Method comparison summary: number of significant sets per test and cohort |\n"
    "| **Table S12** | Effect-size ladder for GSE256068: DEG counts and rates under the four "
    "adjustment models, in the primary comparison and restricted to non-epileptic cortex |",
    "suppl list S12")

# 正文引用 Table S12
rep("Restricted to non-epileptic "
    "cortex the same sequence was 61.3%, 44.3%, 9.5% and 1.0%.",
    "Restricted to non-epileptic "
    "cortex the same sequence was 61.3%, 44.3%, 9.5% and 1.0% (Table S12).",
    "cite S12")

if fail:
    print("!! 失败:", fail); sys.exit(1)
io.open(MS, "w", encoding="utf-8", newline="").write(t)
print("手稿已接入 Table S12")

# ---- 重建 Word ----
r = subprocess.run(["python", os.path.join("scripts",
                                           "101_build_word.py")],
                   capture_output=True, text=True)
print(r.stdout[-600:] if r.returncode == 0 else "构建失败:\n" + r.stderr[:600])

# 把新 docx 拷进包，并同步 md
import shutil
shutil.copy2("03_Manuscript/Manuscript_BIB_BIBformat.docx",
             os.path.join(PKG, "03_Manuscript", "Manuscript_BIBformat.docx"))
shutil.copy2("03_Manuscript/Manuscript_BIB_BIBformat.md",
             os.path.join(PKG, "03_Manuscript", "Manuscript_BIBformat.md"))
print("Word 文件已更新并拷入投稿包")
