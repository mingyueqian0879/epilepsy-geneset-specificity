#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 27: 下载多数据集的表达矩阵
"""
import os, re, json
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "."
DL = os.path.join(BASE, "data", "multi")
os.makedirs(DL, exist_ok=True)
HDR = {"User-Agent": "Mozilla/5.0 (research)"}


def gsm_prefix(acc):
    num = re.match(r"GSE(\d+)$", acc).group(1)
    return "GSE" + (num[:-3] if len(num) > 3 else "0") + "nnn"


FILES = [
    ("GSE134697", "GSE134697_CountMatrix.tsv.gz"),
    ("GSE213488", "GSE213488_grey_matter_count.csv.gz"),
    ("GSE213488", "GSE213488_white_matter_count.csv.gz"),
    ("GSE128300", "GSE128300_RNAstrandedRev_simp_2019.mx.txt.gz"),
    ("GSE272142", "GSE272142_genes.rawcount.anno.xlsx"),
    ("GSE186334", "GSE186334_Cx_cyt_CTLvsW0_countTable_edgeR.txt.gz"),
]


def dl(item):
    acc, fn = item
    url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s"
           % (gsm_prefix(acc), acc, fn))
    dest = os.path.join(DL, fn)
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return (fn, os.path.getsize(dest), "cached")
    try:
        req = urllib.request.Request(url, headers=HDR)
        with urllib.request.urlopen(req, timeout=900) as r, open(dest, "wb") as f:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                f.write(b)
        return (fn, os.path.getsize(dest), "ok")
    except Exception as e:
        return (fn, 0, "ERR " + str(e)[:50])


print("[download] %d files ..." % len(FILES))
with ThreadPoolExecutor(max_workers=4) as ex:
    futs = [ex.submit(dl, f) for f in FILES]
    for f in as_completed(futs):
        fn, sz, st = f.result()
        print("  %-58s %10s  %s" % (fn[:58], "%.1f MB" % (sz / 1e6) if sz else "-", st))

print("\n[files in data/multi]")
for f in sorted(os.listdir(DL)):
    p = os.path.join(DL, f)
    print("  %-58s %8.2f MB" % (f, os.path.getsize(p) / 1e6))
