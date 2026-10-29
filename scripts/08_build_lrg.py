#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 3: 构建 LRG 列表 + ENSG->Symbol 映射
 - 读 Cheng 2023 Table S1 两个 sheet -> 332 LRG
 - 用 Ensembl REST 批量把 LRG symbol 映射成 ENSG
 - 同时用 Gencode v33 / Ensembl 把矩阵的 ENSG 映射成 symbol（用于全基因组背景）
"""
import os, re, json, time, gzip, urllib.request, urllib.parse
import pandas as pd

BASE = "."
OUT  = os.path.join(BASE, "out")
os.makedirs(OUT, exist_ok=True)
HDR  = {"User-Agent": "Mozilla/5.0 (research)", "Content-Type": "application/json"}

# ---------- 1. 读 332 LRG ----------
XL = os.path.join(BASE, "data", "cheng2023_tables", "Table S1 List of 332 lactylation-related genes.xlsx")
sheets = pd.read_excel(XL, sheet_name=None)
print("[Cheng 2023 Table S1]")
lrg = []
for name, df in sheets.items():
    col = df.columns[0]
    vals = [str(v).strip() for v in df[col].dropna().tolist()]
    vals = [v for v in vals if re.fullmatch(r"[A-Za-z0-9\-\.@/]+", v)]
    print(f"  sheet '{name}': {len(vals)} 个")
    print(f"     示例: {', '.join(vals[:12])}")
    lrg += vals

lrg = sorted(set(x.upper().strip() for x in lrg if x))
print(f"\n[合并去重] LRG 总数 = {len(lrg)}")
with open(os.path.join(OUT, "LRG_332.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lrg))
print("[保存] LRG_332.txt")

# ---------- 2. Ensembl REST: symbol -> ENSG ----------
def post_json(url, payload, tries=3):
    data = json.dumps(payload).encode()
    for t in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=HDR, method="POST")
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if t == tries - 1:
                return {"__err__": str(e)}
            time.sleep(2)
    return {}

print("\n[Ensembl REST] LRG symbol -> ENSG ...")
sym2ensg = {}
BATCH = 180
for i in range(0, len(lrg), BATCH):
    chunk = lrg[i:i + BATCH]
    res = post_json("https://rest.ensembl.org/lookup/symbol/homo_sapiens", {"symbols": chunk})
    if "__err__" in res:
        print("  批次失败:", res["__err__"][:70]); continue
    got = 0
    for sym, info in res.items():
        if isinstance(info, dict) and info.get("id"):
            sym2ensg[sym.upper()] = info["id"]        # 可能是多个 -> 取第一个
            got += 1
    print(f"  批次 {i//BATCH+1}: 请求 {len(chunk)} -> 命中 {got}")
    time.sleep(1.0)

print(f"\n[映射结果] {len(sym2ensg)} / {len(lrg)} 个 LRG 拿到 ENSG")
missing = [s for s in lrg if s not in sym2ensg]
print(f"[未映射 {len(missing)}] {', '.join(missing[:35])}")

with open(os.path.join(OUT, "LRG_symbol2ensg.json"), "w", encoding="utf-8") as f:
    json.dump(sym2ensg, f, ensure_ascii=False, indent=1)

# ---------- 3. 保存 LRG 的 ENSG 集合 ----------
lrg_ensg = sorted(set(sym2ensg.values()))
with open(os.path.join(OUT, "LRG_ensg.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lrg_ensg))
print(f"[保存] LRG_ensg.txt  ({len(lrg_ensg)} 个唯一 ENSG)")
