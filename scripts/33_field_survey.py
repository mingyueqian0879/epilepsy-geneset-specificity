#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 33: 领域文献调查 —— 已发表「癫痫基因集生信」论文的报告模式
  检索 PubMed -> 抓摘要 -> 提取报告的 DEG 数 / 基因集 / 是否做特异性检验
"""
import json, os, re, time, urllib.request, urllib.parse
from collections import Counter

OUT = "out"
EUTIL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
HDR = {"User-Agent": "Mozilla/5.0 (research)"}


def esearch(term, retmax=500, **kw):
    q = {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"}
    q.update(kw)
    url = EUTIL + "esearch.fcgi?" + urllib.parse.urlencode(q)
    with urllib.request.urlopen(urllib.request.Request(url, headers=HDR), timeout=120) as r:
        return json.load(r)["esearchresult"]


def efetch_abstracts(pmids):
    """用 efetch 抓摘要"""
    out = {}
    B = 100
    for i in range(0, len(pmids), B):
        chunk = pmids[i:i + B]
        data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(chunk),
                                       "retmode": "text", "rettype": "abstract"}).encode()
        try:
            req = urllib.request.Request(EUTIL + "efetch.fcgi", data=data, headers=HDR)
            with urllib.request.urlopen(req, timeout=180) as r:
                txt = r.read().decode("utf-8", "replace")
        except Exception as e:
            print("   batch fail:", str(e)[:60]); txt = ""
        # 按空行切分记录
        for rec in re.split(r"\n\s*\n(?=\d+\.\s)", txt):
            m = re.match(r"\s*(\d+)\.\s", rec)
            if m:
                idx = int(m.group(1)) - 1
                if idx < len(chunk):
                    out[chunk[idx]] = rec
        time.sleep(0.5)
    return out


# ---------- 检索 ----------
QUERIES = {
 "bioinfo_pipeline": '(epilepsy[Title/Abstract]) AND (bioinformatics[Title/Abstract] OR '
                     '"Gene Expression Omnibus"[Title/Abstract] OR GEO[Title/Abstract] OR '
                     'WGCNA[Title/Abstract] OR "differentially expressed"[Title/Abstract]) '
                     'AND 2020:2026[dp]',
 "geneset_signature": '(epilepsy[Title/Abstract]) AND ("gene signature"[Title/Abstract] OR '
                      '"diagnostic model"[Title/Abstract] OR "hub gene"[Title/Abstract] OR '
                      '"machine learning"[Title/Abstract]) AND ("bioinformatics"[Title/Abstract] '
                      'OR GEO[Title/Abstract] OR WGCNA[Title/Abstract]) AND 2020:2026[dp]',
}
ids = {}
for k, q in QUERIES.items():
    r = esearch(q, retmax=800)
    ids[k] = r["idlist"]
    print("[%s] total=%s fetched=%d" % (k, r.get("count", "?"), len(r["idlist"])))
    time.sleep(0.4)

allids = sorted(set(sum(ids.values(), [])))
print("[unique PMIDs]", len(allids))
abs_map = efetch_abstracts(allids)
print("[abstracts fetched]", len(abs_map))

# ---------- 抽取报告特征 ----------
PAT_DEG = re.compile(r"(\d[\d,]{1,6})\s*(?:differentially expressed genes|DEGs|"
                     r"differentially expressed|genes were differentially)", re.I)
PAT_UP = re.compile(r"(\d[\d,]{1,6})\s*(?:up-?regulated|upregulated)", re.I)
PAT_HUB = re.compile(r"(\d[\d,]{1,6})\s*(?:hub genes|key genes|core genes|hub gene)", re.I)
PAT_TEST = re.compile(r"permutation|negative control|random gene set|specificity test|"
                      r"null distribution|false discovery|type I error", re.I)
PAT_METHOD = re.compile(r"WGCNA|LASSO|SVM-?RFE|random forest|CIBERSORT|ssGSEA|"
                        r"consensus cluster|nomogram|ROC curve", re.I)

recs = []
for pid, txt in abs_map.items():
    title = ""
    m = re.search(r"\n([^\n]{20,300})\n", txt)
    if m:
        title = m.group(1).strip()
    deg = [int(x.replace(",", "")) for x in PAT_DEG.findall(txt)]
    hub = [int(x.replace(",", "")) for x in PAT_HUB.findall(txt)]
    recs.append({
        "pmid": pid,
        "title": title[:200],
        "deg_counts": deg[:5],
        "hub_counts": hub[:5],
        "max_deg": max(deg) if deg else None,
        "has_specificity_test": bool(PAT_TEST.search(txt)),
        "methods": sorted(set(m.group(0).upper() for m in PAT_METHOD.finditer(txt))),
        "len": len(txt),
    })

json.dump(recs, open(os.path.join(OUT, "field_survey_papers.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)

print("\n" + "=" * 76)
print(" 领域调查结果")
print("=" * 76)
n = len(recs)
print("  论文总数                     : %d" % n)
print("  报告了 DEG 数量的论文        : %d (%.1f%%)" %
      (sum(1 for r in recs if r["max_deg"]), 100 * sum(1 for r in recs if r["max_deg"]) / n))
print("  报告了 hub gene 数量的论文   : %d (%.1f%%)" %
      (sum(1 for r in recs if r["hub_counts"]), 100 * sum(1 for r in recs if r["hub_counts"]) / n))
nsp = sum(1 for r in recs if r["has_specificity_test"])
print("  **提到特异性/置换检验的论文  : %d (%.1f%%)**" % (nsp, 100 * nsp / n))
c = Counter()
for r in recs:
    for mm in r["methods"]:
        c[mm] += 1
print("\n  方法学关键词频次:")
for k, v in c.most_common(12):
    print("     %-16s %d (%.1f%%)" % (k, v, 100 * v / n))
degs = sorted(r["max_deg"] for r in recs if r["max_deg"])
if degs:
    import statistics as st
    print("\n  报告的最大 DEG 数: n=%d  中位=%d  均值=%.0f  范围=%d~%d" %
          (len(degs), st.median(degs), st.mean(degs), degs[0], degs[-1]))
    q = [degs[int(len(degs) * p)] for p in (0.25, 0.5, 0.75, 0.9)]
    print("     四分位: Q1=%d  Q2=%d  Q3=%d  P90=%d" % tuple(q))
print("\n[saved] field_survey_papers.json")
