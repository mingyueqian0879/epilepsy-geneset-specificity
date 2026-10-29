#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
追查被反复引用的 332/336 LRG 原始列表来源 + 检查各论文的补充文件
"""
import os, re, json, time, urllib.request, urllib.parse

OUT = "out"
HDR = {"User-Agent": "Mozilla/5.0 (research)"}
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"

def get(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=HDR), timeout=90) as r:
            return r.read().decode("utf-8", "replace")
    except Exception as e:
        return f"__ERR__{e}"

# ---- 1. PMC13502778 : "332 LRGs were obtained from previous research (13)" ----
print("="*76)
print("PMC13502778 : 追查 ref (13) = 332 LRGs 的原始出处")
print("="*76)
xml = get(f"{EPMC}PMC13502778/fullTextXML")
if not xml.startswith("__ERR__"):
    refs = re.findall(r"<ref id=\"([^\"]+)\">(.*?)</ref>", xml, re.S)
    for rid, body in refs:
        n = re.sub(r"<[^>]+>", " ", body)
        n = re.sub(r"\s+", " ", n).strip()
        if re.search(r"\b1[0-9]\b|^13\b", rid) or "13" in rid:
            print(f"  [{rid}] {n[:240]}")
    # 全部列出以便定位
    print("\n  --- 全部参考文献 ---")
    for i, (rid, body) in enumerate(refs, 1):
        n = re.sub(r"<[^>]+>", " ", body)
        n = re.sub(r"\s+", " ", n).strip()
        print(f"  {i:>2}. [{rid}] {n[:190]}")

# ---- 2. PMC13292059 : 336 LRGs from literature [12,13,14,15] ----
print("\n" + "="*76)
print("PMC13292059 : 336 LRGs 来源文献 [12,13,14,15]")
print("="*76)
xml2 = get(f"{EPMC}PMC13292059/fullTextXML")
if not xml2.startswith("__ERR__"):
    refs2 = re.findall(r"<ref id=\"([^\"]+)\">(.*?)</ref>", xml2, re.S)
    for i, (rid, body) in enumerate(refs2, 1):
        if i in (11, 12, 13, 14, 15, 16):
            n = re.sub(r"<[^>]+>", " ", body)
            n = re.sub(r"\s+", " ", n).strip()
            print(f"  {i:>2}. [{rid}] {n[:230]}")

# ---- 3. 检查补充文件 ----
print("\n" + "="*76)
print("检查补充文件")
print("="*76)
for pmcid in ["PMC13502778", "PMC13292059", "PMC13199157", "PMC12219385"]:
    url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}"
    t = get(url)
    tgz = re.findall(r'href="([^"]+)"', t)
    print(f"\n  [{pmcid}]")
    if tgz:
        for h in tgz[:4]:
            print(f"     {h}")
    else:
        print(f"     (无 OA 包)  {t[:120]}")
    time.sleep(0.5)
