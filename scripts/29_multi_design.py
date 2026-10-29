#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 29: 为每个多数据集解析 series matrix -> 样本分组映射
"""
import gzip, io, os, re, json
import urllib.request
from collections import Counter

BASE = "."
OUT = os.path.join(BASE, "out")
HDR = {"User-Agent": "Mozilla/5.0 (research)"}

TARGETS = ["GSE134697", "GSE213488", "GSE128300", "GSE272142", "GSE186334"]


def gsm_prefix(acc):
    num = re.match(r"GSE(\d+)$", acc).group(1)
    return "GSE" + (num[:-3] if len(num) > 3 else "0") + "nnn"


def get_meta(acc):
    url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz"
           % (gsm_prefix(acc), acc, acc))
    req = urllib.request.Request(url, headers=HDR)
    with urllib.request.urlopen(req, timeout=180) as r:
        raw = r.read()
    txt = gzip.decompress(raw).decode("utf-8", "replace")
    tags = {}
    order = []
    for line in txt.split("\n"):
        if line.startswith("!series_matrix_table_begin"):
            break
        if not line.startswith("!Sample_"):
            continue
        p = line.rstrip("\n").split("\t")
        tag = p[0]
        vals = [x.strip().strip('"') for x in p[1:]]
        tags.setdefault(tag, []).append(vals)
        if tag not in order:
            order.append(tag)
    return tags


def build_design(acc, group_keys):
    tags = get_meta(acc)
    n = len(tags["!Sample_geo_accession"][0])
    print("\n" + "=" * 78)
    print("[%s]  n=%d" % (acc, n))
    recs = []
    # 打印所有 characteristics 通道供判断
    for k in tags:
        if "characteristics" in k or k in ("!Sample_title", "!Sample_description",
                                           "!Sample_source_name_ch1"):
            for i, v in enumerate(tags[k]):
                if len(v) != n:
                    continue
                cnt = Counter(v)
                if len(cnt) == 1 and n > 3:
                    continue
                print("   %s [%d] 唯一=%d" % (k, i, len(cnt)))
                for val, c in cnt.most_common(8):
                    print("        %3d  %s" % (c, val[:80]))
    gsm = tags["!Sample_geo_accession"][0]
    title = tags.get("!Sample_title", [[""]] * 1)[0]
    desc = tags.get("!Sample_description", [[""]] * 1)[0]
    recs = {"acc": acc, "n": n, "gsm": gsm, "title": title, "desc": desc}
    return recs, tags


res = {}
for acc in TARGETS:
    try:
        r, tags = build_design(acc, None)
        res[acc] = r
    except Exception as e:
        print("\n[%s] ERROR: %s" % (acc, str(e)[:80]))

json.dump({k: {"acc": v["acc"], "n": v["n"], "gsm": v["gsm"],
               "title": v["title"], "desc": v["desc"]} for k, v in res.items()},
          open(os.path.join(OUT, "multi_design_raw.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\n[saved] multi_design_raw.json")
