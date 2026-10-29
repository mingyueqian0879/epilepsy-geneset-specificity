#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 9: 转向方案可行性检验
  MR/共定位设计不受"细胞组成混杂"影响(工具变量是胚系基因型)。
  关键前提: LRG 基因在脑组织中是否有可用 cis-eQTL 工具变量?
"""
import json, os, time, urllib.request

OUT = "out"
HDR = {"User-Agent": "Mozilla/5.0 (research)", "Accept": "application/json"}

def get(u, tries=3):
    for t in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=HDR), timeout=60) as r:
                return json.load(r)
        except Exception as e:
            if t == tries - 1:
                return {"__err__": str(e)}
            time.sleep(1.5)
    return {}

## 1. GTEx 基因 -> gencodeId 映射
lrg = [l.strip() for l in open(os.path.join(OUT, "LRG_final_symbols.txt"), encoding="utf-8") if l.strip()]
print(f"[LRG] {len(lrg)} 个符号")

# 用 GTEx reference 端点拿 symbol->gencodeId
ref = get("https://gtexportal.org/api/v2/reference/gene?itemsPerPage=1000&geneId=" )
# 改为逐个查询太慢 -> 用 mygene 的 ensembl 结果拼 gencodeId（版本 .x 由 GTEx 自行匹配）
ensg_map = json.load(open(os.path.join(OUT, "LRG_map.json"), encoding="utf-8"))
cand = {}
for sym, v in ensg_map.items():
    cand[sym] = v["ensg"].split(".")[0]
print(f"[映射] {len(cand)} 个 LRG 有 ENSG")

TISSUES = ["Brain_Hippocampus", "Brain_Cortex", "Brain_Frontal_Cortex_BA9"]
TEST = ["GFAP","LDHA","LDHB","EP300","HDAC1","HDAC2","SIRT1","PKM","HIF1A","ALDOA",
        "GAPDH","ENO1","SLC16A1","SLC16A3","CREBBP","SIRT2","SIRT3","H3-3B","TUBA1B","ACTB"]

print(f"\n{'='*74}")
print(" GTEx 脑组织 cis-eQTL 工具变量可得性测试")
print(f"{'='*74}")
ok = 0; tot = 0
for sym in TEST:
    if sym not in cand:
        print(f"  {sym:<9} (无 ENSG)"); continue
    ensg = cand[sym]
    line = f"  {sym:<9} {ensg}"
    found_any = False
    for tis in TISSUES:
        u = (f"https://gtexportal.org/api/v2/association/singleTissueEqtl"
             f"?gencodeId={ensg}&tissueSiteDetailId={tis}&itemsPerPage=3")
        d = get(u)
        n = len(d.get("data", [])) if "__err__" not in d else -1
        line += f"  {tis.replace('Brain_','')}={n}"
        if n and n > 0:
            found_any = True
        time.sleep(0.4)
    tot += 1; ok += 1 if found_any else 0
    print(line + ("   <-- 有工具变量" if found_any else ""))

print(f"\n[小结] 测试 {tot} 个基因, {ok} 个在至少一个脑组织有 cis-eQTL 记录")
print("\n[结论提示] GTEx 脑组织 eQTL 覆盖度直接决定 MR 设计可行性。")
print("           若覆盖率低, 需改用 eQTLGen(全血) 或 BrainMeta(脑 meta-eQTL)。")

## 2. 检查 BrainMeta / 其它脑 eQTL 资源可达性
print(f"\n{'='*74}")
print(" 替代脑 eQTL 资源可达性")
print(f"{'='*74}")
for name, u in [
    ("eQTLGen (全血, 31684样本)", "https://eqtlgen.org/"),
    ("BrainMeta (脑 meta-eQTL)",  "https://hongjiang.info/brainmeta/"),
    ("PsychENCODE",               "https://resource.psychencode.org/"),
    ("OpenGWAS",                  "https://api.opengwas.io/api/status"),
    ("FinnGen",                   "https://r12.finngen.fi/api/"),
]:
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=HDR), timeout=30) as r:
            print(f"  [OK ] {name:<28} HTTP {r.status}")
    except Exception as e:
        print(f"  [??] {name:<28} {str(e)[:60]}")
    time.sleep(0.3)
