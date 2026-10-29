#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量筛查 GEO 癫痫数据集的"分组结构"
方法：只 HTTP Range 抓取 series matrix 前 80KB -> 解压已到达的部分 -> 读 !Sample_characteristics
      （元数据在文件开头，无需下载整个矩阵）
"""
import gzip, io, re, sys, time, json, os
import urllib.request
from collections import Counter

OUT = "out"
os.makedirs(OUT, exist_ok=True)

CANDIDATES = [
    # 常见人类 TLE / 癫痫脑组织候选
    "GSE88992", "GSE201048", "GSE205450", "GSE190452", "GSE202630",
    "GSE245110", "GSE77520", "GSE60772", "GSE100111", "GSE123648",
    "GSE43332", "GSE282281", "GSE16969", "GSE47752", "GSE134697",
    "GSE186334", "GSE143272", "GSE63808", "GSE32550",
    "GSE110513", "GSE115432", "GSE103156", "GSE137217",
]

HDR = {"User-Agent": "Mozilla/5.0 (research; GEO metadata screening)"}

def gsm_prefix(acc):
    """GSE12345 -> GSE12nnn  (GEO FTP 分目录规则)"""
    m = re.match(r"(GSE)(\d+)$", acc)
    num = m.group(2)
    stem = num[:-3] if len(num) > 3 else "0"
    return f"GSE{stem}nnn"

def peek_metadata(acc, nbytes=80000):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gsm_prefix(acc)}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    req = urllib.request.Request(url, headers={**HDR, "Range": f"bytes=0-{nbytes}"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            raw = r.read()
    except Exception as e:
        return {"acc": acc, "error": f"download: {e}"}
    # 解压已收到的部分（尾部可能不完整 -> 忽略异常）
    try:
        txt = gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode("utf-8", "replace")
    except Exception:
        try:
            d = gzip.GzipFile(fileobj=io.BytesIO(raw))
            chunks = []
            while True:
                b = d.read(65536)
                if not b: break
                chunks.append(b)
            txt = b"".join(chunks).decode("utf-8", "replace")
        except Exception:
            try:
                txt = gzip.decompress(raw).decode("utf-8", "replace")
            except Exception:
                txt = raw.decode("utf-8", "replace")

    title = ""
    chars = []
    n_title = 0
    for line in txt.split("\n"):
        if line.startswith("!Series_title"):
            title = line.split("\t", 1)[-1].strip().strip('"')[:110]
        elif line.startswith("!Sample_title"):
            n_title = len(line.rstrip().split("\t")) - 1
        elif line.startswith("!Sample_characteristics"):
            chars.append([c.strip().strip('"') for c in line.rstrip().split("\t")[1:]])
        if line.startswith("!series_matrix_table_begin"):
            break
    return {"acc": acc, "title": title, "n": n_title, "chars": chars}

results = []
for acc in CANDIDATES:
    r = peek_metadata(acc)
    results.append(r)
    if "error" in r:
        print(f"[ERR ] {acc:<12} {r['error'][:70]}")
    else:
        print(f"\n[ OK ] {acc:<12} n={r['n']}")
        print(f"       {r['title']}")
        for i, c in enumerate(r["chars"]):
            cnt = Counter(c)
            if len(cnt) == 1 and len(c) == r["n"]:
                continue  # 常量字段跳过
            print(f"       char[{i}] 唯一={len(cnt)}")
            for v, k in cnt.most_common(10):
                print(f"           {k:>4}  {v[:85]}")
    time.sleep(0.6)

with open(os.path.join(OUT, "geo_screen.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=1)
print(f"\n[保存] {os.path.join(OUT,'geo_screen.json')}")
