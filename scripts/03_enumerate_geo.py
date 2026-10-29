#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
正规做法：先枚举 GEO 全部人类癫痫表达数据集（含样本数），再定向抓表头确认分组
"""
import gzip, io, re, os, json, time, urllib.request, urllib.parse
from collections import Counter

OUT = "out"
os.makedirs(OUT, exist_ok=True)
EUTIL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HDR = {"User-Agent": "Mozilla/5.0 (research)"}

def esearch(db, term, retmax=500, **kw):
    q = {"db": db, "term": term, "retmax": retmax, "retmode": "json"}
    q.update(kw)
    url = EUTIL + "esearch.fcgi?" + urllib.parse.urlencode(q)
    with urllib.request.urlopen(urllib.request.Request(url, headers=HDR), timeout=90) as r:
        return json.load(r)["esearchresult"]

def esummary(db, ids):
    url = EUTIL + "esummary.fcgi"
    data = urllib.parse.urlencode({"db": db, "id": ",".join(map(str, ids)), "retmode": "json"}).encode()
    with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=HDR), timeout=120) as r:
        return json.load(r)["result"]

# ---------- 1. 枚举人类癫痫 GEO 系列 ----------
print("=" * 70)
print("步骤1: 枚举 GEO 人类癫痫表达数据集")
print("=" * 70)
term = ('(epilepsy[All Fields] OR "temporal lobe epilepsy"[All Fields] OR "seizure"[All Fields]) '
        'AND "Homo sapiens"[Organism] AND gse[Entry Type]')
res = esearch("gds", term, retmax=500)
total = int(res["count"]); ids = res["idlist"]
print(f"总命中 = {total}，本轮取 {len(ids)}")
time.sleep(0.4)

rows = []
for i in range(0, len(ids), 100):
    batch = ids[i:i + 100]
    s = esummary("gds", batch)
    for uid in batch:
        d = s.get(uid)
        if not d or d.get("entrytype") != "GSE":
            continue
        try:
            n = int(d.get("n_samples") or 0)
        except Exception:
            n = 0
        rows.append({
            "acc": d.get("accession", ""),
            "n": n,
            "year": (d.get("pdat") or "").split("/")[0],
            "title": (d.get("title") or "").replace("\n", " "),
            "summary": (d.get("summary") or "").replace("\n", " ")[:300],
            "gpl": d.get("gpl", ""),
            "gdstype": d.get("gdstype", ""),
        })
    time.sleep(0.4)

print(f"解析到 GSE 记录 = {len(rows)}")
with open(os.path.join(OUT, "geo_all_epilepsy.json"), "w", encoding="utf-8") as f:
    json.dump(rows, f, ensure_ascii=False, indent=1)

# ---------- 2. 筛选候选：脑组织 + 样本数够 + 含 control 线索 ----------
print("\n" + "=" * 70)
print("步骤2: 筛选候选（n>=15 且标题/摘要提示脑组织或对照）")
print("=" * 70)
BRAIN = re.compile(r"hippocamp|cortex|cortical|brain|temporal lobe|TLE|epilep", re.I)
CTRL  = re.compile(r"control|postmortem|post-mortem|non-epilep|nonepilep|normal|healthy", re.I)

cand = []
for r in rows:
    if r["n"] < 15:
        continue
    blob = r["title"] + " " + r["summary"]
    if not BRAIN.search(blob):
        continue
    r["has_ctrl_hint"] = bool(CTRL.search(blob))
    cand.append(r)

cand.sort(key=lambda x: (not x["has_ctrl_hint"], -x["n"]))
print(f"候选数 = {len(cand)}\n")
for r in cand[:45]:
    flag = "CTRL?" if r["has_ctrl_hint"] else "     "
    print(f"  {flag} {r['acc']:<12} n={r['n']:<4} {r['year']}  {r['title'][:88]}")

with open(os.path.join(OUT, "geo_candidates.json"), "w", encoding="utf-8") as f:
    json.dump(cand, f, ensure_ascii=False, indent=1)
print(f"\n[保存] geo_candidates.json ({len(cand)} 条)")
