#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 129: 对已交付图件做「文字包围盒」确定性重叠检测。

方法：tesseract TSV 取词级 bbox（置信度 >= 40）。
关键控制：
  * 只在**同一次 psm 运行内**两两比较，避免两次识别同一文字造成假重叠；
  * 要求两个词文本不同（同词重复框视为同一文字）；
  * 交叠率 inter / min(area) > 0.30 判为真重叠。
另检贴边（可能被裁切）。
这一层只看"文字压文字"；"文字压图形"由目视 + 定向裁切复核。
"""
import os
import sys
import subprocess
import itertools
import json

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

D = "figures"
FIGS = ["Fig1_field_survey.png", "Fig3_Z_distribution.png",
        "Fig4_hot_vs_positive.png", "Fig5_label_permutation.png",
        "Fig6_cell_composition.png", "Fig7_reporting_standard.png",
        "Fig8_reproducibility.png", "Fig9_spikein_validation.png",
        "graphical_abstract.png"]


def ocr_words(path, psm):
    out = subprocess.run(["tesseract", path, "stdout", "-l", "eng",
                          "--psm", str(psm), "-c", "tessedit_create_tsv=1"],
                         capture_output=True, text=True, encoding="utf-8",
                         errors="replace", timeout=240)
    rows = []
    for ln in out.stdout.splitlines()[1:]:
        f = ln.split("\t")
        if len(f) < 12:
            continue
        try:
            conf = float(f[10])
        except ValueError:
            continue
        txt = f[11].strip()
        if not txt or conf < 40:
            continue
        left, top, w, h = int(f[6]), int(f[7]), int(f[8]), int(f[9])
        if w <= 0 or h <= 0:
            continue
        rows.append(dict(left=left, top=top, w=w, h=h, conf=conf, text=txt,
                         line=(f[2], f[3], f[4])))
    return rows


def inter_ratio(a, b):
    ax0, ay0, ax1, ay1 = a["left"], a["top"], a["left"] + a["w"], a["top"] + a["h"]
    bx0, by0, bx1, by1 = b["left"], b["top"], b["left"] + b["w"], b["top"] + b["h"]
    ix = max(0, min(ax1, bx1) - max(ax0, bx0))
    iy = max(0, min(ay1, by1) - max(ay0, by0))
    if ix == 0 or iy == 0:
        return 0.0
    return (ix * iy) / float(min(a["w"] * a["h"], b["w"] * b["h"]))


def dedup(rows):
    """同一文字在同一位置附近的重复框只保留一个。"""
    keep = []
    for r in sorted(rows, key=lambda x: -x["conf"]):
        dup = False
        for k in keep:
            if k["text"] == r["text"] and inter_ratio(k, r) > 0.5:
                dup = True
                break
        if not dup:
            keep.append(r)
    return keep


def safe(s):
    return "".join(ch if 32 <= ord(ch) < 127 else "?" for ch in s)


def audit(name):
    path = os.path.join(D, name)
    print("=" * 78)
    print(name)
    print("=" * 78)
    hits_total, edge_total, words_total = [], [], 0
    for psm in (3, 11):
        rows = dedup(ocr_words(path, psm))
        words_total = max(words_total, len(rows))
        hits = []
        for a, b in itertools.combinations(rows, 2):
            if a["text"] == b["text"]:
                continue
            r = inter_ratio(a, b)
            if r > 0.30:
                hits.append((round(r, 2), a["text"], b["text"],
                             (a["left"], a["top"], a["w"], a["h"]),
                             (b["left"], b["top"], b["w"], b["h"])))
        hits.sort(reverse=True)
        print("  psm=%-2d words=%d  overlapping pairs=%d" % (psm, len(rows), len(hits)))
        for h in hits[:10]:
            print("     %.2f  %-20s <-> %-20s  A%s B%s"
                  % (h[0], safe(h[1]), safe(h[2]), h[3], h[4]))
        hits_total.extend(hits)
        from PIL import Image
        W, H = Image.open(path).size
        edge = [r for r in rows
                if r["left"] <= 1 or r["top"] <= 1
                or r["left"] + r["w"] >= W - 1 or r["top"] + r["h"] >= H - 1]
        edge_total.extend(edge)
        print("     boxes touching border: %d%s"
              % (len(edge), "".join(" [%s]" % safe(e["text"]) for e in edge[:5])))
    print()
    return {"figure": name, "words": words_total,
            "overlaps": [h[:3] for h in hits_total], "n_edge": len(edge_total)}


if __name__ == "__main__":
    report = [audit(f) for f in FIGS]
    print("#" * 78)
    print("SUMMARY  (overlap = word boxes from the same OCR pass intersecting >30%)")
    print("#" * 78)
    for r in report:
        n = len(r["overlaps"])
        print("  %-34s words=%-4d overlaps=%-3d edge=%-3d%s"
              % (r["figure"], r["words"], n, r["n_edge"],
                 "   <== CHECK" if n else ""))
    with open("out/_layout_ocr_audit.json", "w",
              encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
