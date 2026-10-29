#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 25: Batch screen for human epilepsy datasets that support case-control design.
  Fetch series-matrix header via HTTP Range, auto-detect grouping fields.
"""
import gzip, io, os, re, json
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter

BASE = "."
OUT = os.path.join(BASE, "out")
HDR = {"User-Agent": "Mozilla/5.0 (research)"}

rows = json.load(open(os.path.join(OUT, "geo_all_epilepsy.json"), encoding="utf-8"))
print("[input] enumerated records =", len(rows))

BRAIN = re.compile(r"hippocamp|cortex|cortical|brain|temporal lobe|TLE|epilep|seizure", re.I)
cands = {r["acc"]: r for r in rows
         if r["n"] >= 12 and BRAIN.search(r["title"] + " " + r["summary"])}
print("[candidates] n>=12 and brain/epilepsy related =", len(cands))


def gsm_prefix(acc):
    num = re.match(r"GSE(\d+)$", acc).group(1)
    return "GSE" + (num[:-3] if len(num) > 3 else "0") + "nnn"


def peek(acc, nbytes=800000):
    url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz"
           % (gsm_prefix(acc), acc, acc))
    req = urllib.request.Request(url, headers=dict(HDR, Range="bytes=0-%d" % nbytes))
    try:
        with urllib.request.urlopen(req, timeout=90) as r:
            raw = r.read()
    except Exception as e:
        return {"acc": acc, "err": str(e)[:60], "n": 0}

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
        try:
            txt = gzip.decompress(raw).decode("utf-8", "replace")
        except Exception:
            txt = raw.decode("utf-8", "replace")

    title = ""
    chars = []
    ns = 0
    plat = ""
    for line in txt.split("\n"):
        if line.startswith("!Series_title"):
            title = line.split("\t", 1)[-1].strip().strip('"')[:150]
        elif line.startswith("!Series_platform_id"):
            plat = line.split("\t", 1)[-1].strip().strip('"')
        elif line.startswith("!Sample_title"):
            ns = len(line.rstrip().split("\t")) - 1
        elif line.startswith("!Sample_characteristics"):
            chars.append([c.strip().strip('"') for c in line.rstrip().split("\t")[1:]])
        elif line.startswith("!series_matrix_table_begin"):
            break

    group_field = None
    groups = {}
    for i, c in enumerate(chars):
        if ns == 0 or len(c) != ns:
            continue
        cnt = Counter(c)
        if len(cnt) < 2 or len(cnt) > 30:
            continue
        blob = " ".join(cnt.keys()).lower()
        has_case = bool(re.search(r"epilep|patient|tle|seizure|case|hs\b|fcd|tsc|disease", blob))
        has_ctrl = bool(re.search(r"control|normal|healthy|non-epilep|nonepilep|without epilepsy", blob))
        if has_case and has_ctrl:
            group_field = "char[%d]" % i
            groups = dict(cnt)
            break
    return {"acc": acc, "title": title, "n": ns, "plat": plat,
            "group_field": group_field, "groups": groups, "n_char": len(chars)}


accs = sorted(cands.keys())
print("[fetch headers]", len(accs), "datasets ...")
res = []
with ThreadPoolExecutor(max_workers=6) as ex:
    futs = {ex.submit(peek, a): a for a in accs}
    for k, f in enumerate(as_completed(futs), 1):
        res.append(f.result())
        if k % 25 == 0:
            print("   ... %d/%d" % (k, len(accs)))
print("[done]", len(res))

ok = [r for r in res if r.get("group_field")]
ok.sort(key=lambda x: -x["n"])
print("=" * 80)
print("FOUND %d datasets with case+control grouping" % len(ok))
print("=" * 80)
for r in ok:
    gs = "; ".join("%s x%d" % (v, c) for v, c in
                   sorted(r["groups"].items(), key=lambda t: -t[1])[:6])
    print("\n[%s]  n=%d  %s" % (r["acc"], r["n"], r["plat"]))
    print("   " + r["title"][:110])
    print("   field %s : %s" % (r["group_field"], gs))

json.dump(res, open(os.path.join(OUT, "dataset_screen2.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
json.dump(ok, open(os.path.join(OUT, "dataset_casectrl.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\n[saved] dataset_screen2.json / dataset_casectrl.json")