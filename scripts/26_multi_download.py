#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 26: 下载多数据集 —— 检查每个的表达数据获取方式（矩阵内 vs 补充文件）
"""
import gzip, io, os, re, json
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "."
OUT = os.path.join(BASE, "out")
DL = os.path.join(BASE, "data", "multi")
os.makedirs(DL, exist_ok=True)
HDR = {"User-Agent": "Mozilla/5.0 (research)"}

TARGETS = ["GSE186334", "GSE168375", "GSE213488", "GSE272142",
           "GSE128300", "GSE245228", "GSE134697", "GSE256068"]


def gsm_prefix(acc):
    num = re.match(r"GSE(\d+)$", acc).group(1)
    return "GSE" + (num[:-3] if len(num) > 3 else "0") + "nnn"


def fetch(url, nbytes=None, timeout=300):
    hdr = dict(HDR)
    if nbytes:
        hdr["Range"] = "bytes=0-%d" % nbytes
    req = urllib.request.Request(url, headers=hdr)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def probe(acc):
    base = ("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/" % (gsm_prefix(acc), acc))
    info = {"acc": acc}
    # 1) series matrix
    sm = base + "matrix/%s_series_matrix.txt.gz" % acc
    try:
        raw = fetch(sm, nbytes=600000)
        try:
            d = gzip.GzipFile(fileobj=io.BytesIO(raw))
            buf = []
            while True:
                b = d.read(131072)
                if not b:
                    break
                buf.append(b)
            txt = b"".join(buf).decode("utf-8", "replace")
        except Exception:
            txt = ""
        # 是否有表达表
        has_table = False
        nrow = 0
        in_tbl = False
        for line in txt.split("\n"):
            if line.startswith("!series_matrix_table_begin"):
                in_tbl = True
                continue
            if in_tbl:
                if line.startswith("!series_matrix_table_end"):
                    break
                nrow += 1
        has_table = nrow > 100
        info["matrix_rows"] = nrow
        # 平台
        m = re.search(r'!Series_platform_id\t"?([A-Z0-9]+)', txt)
        info["plat"] = m.group(1) if m else ""
        # 补充文件
        supp = re.findall(r'!Series_supplementary_file\t"?([^\t\n"]+)', txt)
        if not supp:
            supp = re.findall(r'!Series_supplementary_file_\d+\t"?([^\t\n"]+)', txt)
        info["series_suppl"] = [s.strip() for s in supp][:5]
    except Exception as e:
        info["sm_err"] = str(e)[:60]

    # 2) 目录下的 supplementary 列表
    try:
        page = fetch(base + "suppl/", timeout=60).decode("utf-8", "replace")
        files = re.findall(r'href="([^"?][^"]*\.(?:gz|txt|csv|tsv|xlsx|tar|zip))"', page)
        info["suppl_files"] = sorted(set(files))[:10]
    except Exception as e:
        info["suppl_err"] = str(e)[:50]
    return info


print("[probe] %d datasets ..." % len(TARGETS))
res = []
with ThreadPoolExecutor(max_workers=4) as ex:
    futs = {ex.submit(probe, a): a for a in TARGETS}
    for f in as_completed(futs):
        res.append(f.result())

for r in sorted(res, key=lambda x: x["acc"]):
    print("\n" + "=" * 74)
    print("[%s]  platform=%s  matrix_rows=%s" % (r["acc"], r.get("plat", "?"), r.get("matrix_rows", "?")))
    if r.get("sm_err"):
        print("   series_matrix error:", r["sm_err"])
    if r.get("suppl_files"):
        for f in r["suppl_files"]:
            print("   suppl:", f)
    else:
        print("   suppl: (none)")

json.dump(res, open(os.path.join(OUT, "multi_probe.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\n[saved] multi_probe.json")
