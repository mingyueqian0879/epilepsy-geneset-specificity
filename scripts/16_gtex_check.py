#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""GTEx API 正确用法校正测试"""
import json, time, urllib.request

HDR = {"User-Agent": "Mozilla/5.0 (research)", "Accept": "application/json"}

def get(u):
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=HDR), timeout=60) as r:
            return json.load(r)
    except Exception as e:
        return {"__err__": str(e)}

print("=== 1) reference/gene 拿正确 gencodeId ===")
ids = {}
for sym in ["GFAP", "GAPDH", "LDHA", "EP300", "SIRT1"]:
    d = get(f"https://gtexportal.org/api/v2/reference/gene?geneId={sym}")
    if "__err__" in d:
        print(f"  {sym}: ERR {d['__err__'][:80]}")
    else:
        dat = d.get("data", [])
        if dat:
            x = dat[0]
            ids[sym] = x.get("gencodeId")
            print(f"  {sym:<7} gencodeId={x.get('gencodeId')}  ensembl={x.get('ensemblGeneId')}")
        else:
            print(f"  {sym}: 空, keys={list(d.keys())}")
    time.sleep(0.4)

print("\n=== 2) 用正确 gencodeId 查脑 eQTL ===")
TIS = ["Brain_Cortex", "Brain_Hippocampus", "Brain_Frontal_Cortex_BA9"]
for sym, gid in ids.items():
    if not gid:
        continue
    row = f"  {sym:<7}"
    for tis in TIS:
        d = get(f"https://gtexportal.org/api/v2/association/singleTissueEqtl"
                f"?gencodeId={gid}&tissueSiteDetailId={tis}&itemsPerPage=5")
        n = len(d.get("data", [])) if "__err__" not in d else -1
        row += f"  {tis.replace('Brain_','')}={n}"
        if n and n > 0 and tis == "Brain_Cortex":
            r0 = d["data"][0]
            row += f"(p={r0.get('pval'):.1e})"
        time.sleep(0.35)
    print(row)

print("\n=== 3) 结论 ===")
print("  若上面有非0记录 -> GTEx 脑 eQTL API 可用, MR 设计可行")
print("  若仍全0        -> 需下载 GTEx v8 全量 eQTL 文件 (不等API)")
