#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 3b: 稳健构建 LRG 集 + 全矩阵 ENSG<->Symbol 双向映射 (mygene.info, 别名感知)
"""
import os, re, json, time, urllib.request
import pandas as pd

BASE = "."
OUT  = os.path.join(BASE, "out")
HDR  = {"User-Agent": "Mozilla/5.0 (research)", "Content-Type": "application/json"}

def mygene(q_list, scopes="symbol,alias", fields="ensembl.gene,symbol", species="human", tries=3):
    """批量查询 mygene.info；返回 {query_term: {ensg, symbol}}"""
    url = "https://mygene.info/v3/query"
    out = {}
    B = 400
    for i in range(0, len(q_list), B):
        chunk = q_list[i:i + B]
        payload = {"q": chunk, "scopes": scopes, "fields": fields,
                   "species": species, "size": 1}
        for t in range(tries):
            try:
                req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                             headers=HDR, method="POST")
                with urllib.request.urlopen(req, timeout=120) as r:
                    res = json.load(r)
                break
            except Exception as e:
                if t == tries - 1:
                    print(f"    批次 {i//B+1} 失败: {e}"); res = []
                time.sleep(2)
        if isinstance(res, dict):   # 单条结果时 mygene 返回 dict
            res = [res]
        for item in res or []:
            if not isinstance(item, dict):
                continue
            q = str(item.get("query", "")).upper()
            ens = item.get("ensembl")
            ensg = None
            if isinstance(ens, dict):
                ensg = ens.get("gene")
            elif isinstance(ens, list) and ens:
                ensg = ens[0].get("gene") if isinstance(ens[0], dict) else None
            if q and ensg:
                out[q] = {"ensg": ensg, "symbol": item.get("symbol", "")}
        print(f"    批次 {i//B+1}: 请求 {len(chunk)}  累计命中 {len(out)}")
        time.sleep(0.6)
    return out

# ---------- 1. 读完整 332 LRG ----------
XL = os.path.join(BASE, "data", "cheng2023_tables", "Table S1 List of 332 lactylation-related genes.xlsx")
sheets = pd.read_excel(XL, sheet_name=None)
lrg = []
for name, df in sheets.items():
    for col in df.columns:
        vals = [str(v).strip() for v in df[col].dropna().tolist()]
        vals = [v for v in vals if re.fullmatch(r"[A-Za-z0-9\-\.@/]+", v)]
        print(f"  sheet '{name}' 列 '{col}': {len(vals)}")
        lrg += vals
lrg = sorted(set(x.upper() for x in lrg))
print(f"\n[LRG 合并去重] = {len(lrg)}")
with open(os.path.join(OUT, "LRG_final_symbols.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lrg))

# ---------- 2. LRG symbol -> ENSG ----------
print("\n[mygene] LRG -> ENSG ...")
lrg_map = mygene(lrg)
print(f"  LRG 命中 = {len(lrg_map)} / {len(lrg)}")
miss = [s for s in lrg if s not in lrg_map]
print(f"  未命中({len(miss)}): {', '.join(miss)}")

lrg_ensg = sorted({v["ensg"].split(".")[0] for v in lrg_map.values()})
with open(os.path.join(OUT, "LRG_final_ensg.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lrg_ensg))
with open(os.path.join(OUT, "LRG_map.json"), "w", encoding="utf-8") as f:
    json.dump(lrg_map, f, ensure_ascii=False, indent=1)
print(f"[保存] LRG_final_symbols.txt ({len(lrg)})  LRG_final_ensg.txt ({len(lrg_ensg)})")

# ---------- 3. 全矩阵 ENSG -> symbol ----------
expr_genes = pd.read_pickle(os.path.join(OUT, "expr_log2tpm.pkl")).index.tolist()
ensg_list = [g.split(".")[0] for g in expr_genes]
print(f"\n[mygene] 全矩阵 {len(ensg_list)} 个 ENSG -> symbol ...")
full = mygene(ensg_list, scopes="ensembl.gene", fields="symbol,ensembl.gene")
print(f"  全矩阵命中 = {len(full)}")
with open(os.path.join(OUT, "ensg2symbol.json"), "w", encoding="utf-8") as f:
    json.dump(full, f, ensure_ascii=False)
print("[保存] ensg2symbol.json")
