#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
从 Europe PMC 开放获取的乳酸化生信论文中提取其实际使用的 LRG 基因列表
"""
import os, re, json, sys, time, urllib.request, urllib.parse

OUT = "out"
os.makedirs(OUT, exist_ok=True)
HDR = {"User-Agent": "Mozilla/5.0 (research)"}
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"

def fulltext_xml(pmcid):
    url = f"{EPMC}{pmcid}/fullTextXML"
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=HDR), timeout=90) as r:
            return r.read().decode("utf-8", "replace")
    except Exception as e:
        return f"__ERR__{e}"

def search(q, n=25):
    url = EPMC + "search?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "pageSize": n, "resultType": "core"})
    with urllib.request.urlopen(urllib.request.Request(url, headers=HDR), timeout=90) as r:
        return json.load(r)

# 目标：明确用了 GeneCards / 提供基因列表的乳酸化生信论文
q = ('TITLE_ABS:"lactylation" AND (TITLE_ABS:"bioinformatics" OR TITLE_ABS:"WGCNA" '
     'OR TITLE_ABS:"gene signature" OR TITLE_ABS:"diagnostic")')
res = search(q, 40)
items = [x for x in res["resultList"]["result"] if x.get("isOpenAccess") == "Y" and x.get("pmcid")]
print(f"开放获取候选 = {len(items)}\n")

GENE_RE = re.compile(r"\b[A-Z][A-Z0-9]{1,9}\b")
KEY = re.compile(r"lactylation[- ]related gene|lactylation[- ]associated gene|LRG|GeneCards|"
                 r"lactylation regulator|gene set", re.I)

found = []
for it in items[:12]:
    pmcid = it["pmcid"]
    title = it["title"][:95]
    xml = fulltext_xml(pmcid)
    if xml.startswith("__ERR__"):
        print(f"[ERR] {pmcid} {xml[:60]}")
        time.sleep(0.5)
        continue
    print(f"\n{'='*74}\n[{pmcid}] {title}")
    # 找提到基因集来源的句子
    plain = re.sub(r"<[^>]+>", " ", xml)
    plain = re.sub(r"\s+", " ", plain)
    for m in KEY.finditer(plain):
        s = max(0, m.start() - 200); e = min(len(plain), m.end() + 260)
        snippet = plain[s:e]
        if re.search(r"GeneCards|MSigDB|obtained|retrieved|downloaded|collected|list of", snippet, re.I):
            print("   >>", snippet.strip()[:430])
            break
    # 统计正文中出现的基因 symbol 数（粗筛）
    supp = xml.lower().find("supplementary")
    print(f"   正文长度={len(plain)}  含supplementary={supp>0}")
    found.append({"pmcid": pmcid, "title": title})
    time.sleep(0.6)

with open(os.path.join(OUT, "lrg_source_papers.json"), "w", encoding="utf-8") as f:
    json.dump(found, f, ensure_ascii=False, indent=1)
print(f"\n[保存] lrg_source_papers.json")
