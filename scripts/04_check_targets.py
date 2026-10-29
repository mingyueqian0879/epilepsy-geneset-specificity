#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
定向核查重点候选数据集的"分组结构"（加大抓取范围到 600KB）
"""
import gzip, io, re, os, json, time, urllib.request
from collections import Counter

OUT = "out"
HDR = {"User-Agent": "Mozilla/5.0 (research)"}

TARGETS = [
    "GSE256068",   # n=162 drug-resistant epilepsies  <- 重点
    "GSE190453",   # n=79  epilepsy astrocytes (bulk)
    "GSE163296",   # n=22  hippocampal CA3
    "GSE121010",   # n=24  Rasmussen vs control
    "GSE310367",   # n=49  Rasmussen encephalitis
    "GSE99455",    # n=24  hippocampal miRNA
    "GSE168916",   # n=33  epilepsy DNA methylation
    "GSE186334",   # n=68  mTLE + postmortem controls (已确认有对照)
    "GSE134697",   # n=36  mTLE (2 对照)
    "GSE119684",   # n=45  epilepsy methylome
    "GSE156374",   # n=96  polymicrogyria
]

def gsm_prefix(acc):
    num = re.match(r"GSE(\d+)$", acc).group(1)
    return f"GSE{num[:-3] if len(num)>3 else 0}nnn"

def peek(acc, nbytes=600000):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gsm_prefix(acc)}/{acc}/matrix/{acc}_series_matrix.txt.gz"
    req = urllib.request.Request(url, headers={**HDR, "Range": f"bytes=0-{nbytes}"})
    try:
        with urllib.request.urlopen(req, timeout=90) as r:
            raw = r.read()
    except Exception as e:
        return {"acc": acc, "err": f"{e}"}
    txt = ""
    try:
        d = gzip.GzipFile(fileobj=io.BytesIO(raw))
        buf = []
        while True:
            b = d.read(131072)
            if not b: break
            buf.append(b)
        txt = b"".join(buf).decode("utf-8", "replace")
    except Exception:
        try:
            txt = gzip.decompress(raw).decode("utf-8", "replace")
        except Exception as e:
            return {"acc": acc, "err": f"gunzip: {e}"}
    title = ""; chars = []; npr = 0; plats = set(); nsamp = 0
    for line in txt.split("\n"):
        if line.startswith("!Series_title"):        title = line.split("\t",1)[-1].strip().strip('"')[:130]
        elif line.startswith("!Series_platform_id"): plats.add(line.split("\t",1)[-1].strip().strip('"'))
        elif line.startswith("!Sample_title"):       nsamp = len(line.rstrip().split("\t"))-1
        elif line.startswith("!Sample_characteristics"):
            chars.append([c.strip().strip('"') for c in line.rstrip().split("\t")[1:]])
        elif line.startswith("!series_matrix_table_begin"):
            break
    return {"acc": acc, "title": title, "n": nsamp, "plats": sorted(plats), "chars": chars}

out = []
for acc in TARGETS:
    r = peek(acc)
    out.append(r)
    print("\n" + "="*72)
    if "err" in r:
        print(f"[ERR] {acc}: {r['err'][:90]}")
    else:
        print(f"[{acc}]  n={r['n']}  platform={r.get('plats')}")
        print(f"  {r['title']}")
        if r["n"] == 0:
            print("  (!) 元数据超出抓取范围")
        for i, c in enumerate(r["chars"]):
            cnt = Counter(c)
            if len(cnt) == 1 and len(c) == r["n"]:
                print(f"  char[{i}] 常量: {list(cnt)[0][:70]}")
                continue
            print(f"  char[{i}] 唯一={len(cnt)}")
            for v, k in cnt.most_common(12):
                print(f"      {k:>4}  {v[:82]}")
    time.sleep(0.8)

with open(os.path.join(OUT, "geo_targets.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print(f"\n[保存] geo_targets.json")
