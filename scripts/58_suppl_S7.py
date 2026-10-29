#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 58: 生成补充表 S7 (逐基因效力) + 补充方法, 并重建投稿包"""
import os, json
import pandas as pd

PKG = "."
SUP = os.path.join(PKG, "06_Supplementary")
OUT = "out"
os.makedirs(SUP, exist_ok=True)

D = json.load(open(os.path.join(OUT, "P0_revision_numbers.json"), encoding="utf-8"))

# ---------- Table S7: 逐基因效力 ----------
s7 = pd.DataFrame(D["power"])
s7.columns = ["Cohort", "Cases", "Controls", "Power (%) d=0.5",
              "Power (%) d=0.8", "Power (%) d=1.0"]
s7 = s7.round(1)
s7.to_csv(os.path.join(SUP, "TableS7_statistical_power.csv"), index=False)

# ---------- Table S8: 最小可检出富集 ----------
s8 = pd.DataFrame(D["mde"])
if len(s8):
    s8 = s8[["dataset", "m", "deg_rate", "z_need", "exp_k100", "rel_excess", "extra_deg"]]
    s8.columns = ["Cohort", "Gene sets tested", "DEG rate (%)",
                  "Nominal P required (Z)", "Expected DEGs in 100-gene set",
                  "Minimum detectable relative excess (%)",
                  "Extra DEGs required"]
    s8 = s8.round(3)
    s8.to_csv(os.path.join(SUP, "TableS8_minimum_detectable_enrichment.csv"), index=False)

# ---------- Table S9: 校准探针排名 ----------
rows = []
for ds, items in D["positive_controls"].items():
    for it in items:
        rows.append({"Cohort": ds, "Calibration probe": it["set"],
                     "Genes detected": it["n"], "Observed DEGs": it["obs"],
                     "Expected DEGs": round(it["exp"], 2), "Z": round(it["Z"], 2),
                     "Rank": it["rank"], "Sets tested": it["m"],
                     "Rank percentile": round(it["pct"], 1)})
s9 = pd.DataFrame(rows)
if len(s9):
    s9.to_csv(os.path.join(SUP, "TableS9_calibration_probes.csv"), index=False)

# ---------- 补充方法: 完整检索式 ----------
sm = """# Supplementary Methods

## S1. Full PubMed search strategies for the literature survey

Queried 22 September 2026.

**Strategy (i)**

```
(epilepsy[Title/Abstract]) AND (bioinformatics[Title/Abstract] OR
"Gene Expression Omnibus"[Title/Abstract] OR GEO[Title/Abstract] OR
WGCNA[Title/Abstract] OR "differentially expressed"[Title/Abstract])
AND 2020:2026[dp]
```

**Strategy (ii)**

```
(epilepsy[Title/Abstract]) AND ("gene signature"[Title/Abstract] OR
"diagnostic model"[Title/Abstract] OR "hub gene"[Title/Abstract] OR
"machine learning"[Title/Abstract]) AND (bioinformatics[Title/Abstract] OR
GEO[Title/Abstract] OR WGCNA[Title/Abstract]) AND 2020:2026[dp]
```

Records were deduplicated to 582 unique PMIDs; abstracts were retrieved for 576 via
NCBI EFetch. Screening was performed on abstracts only, as stated in the main text,
and therefore provides a lower bound on the true prevalence of specificity testing.

## S2. Computation of the minimum detectable enrichment (MDE)

For a cohort with `m` gene sets tested and a genome-wide DEG proportion `p`, a gene
set of `k` detected genes has an expression-matched expectation of `k·p` differentially
expressed genes. Under the most favourable case — the gene set attaining the smallest
P value, i.e. rank 1 in the Benjamini–Hochberg ordering — the nominal significance
threshold is `0.05/m`. The corresponding Z requirement is

```
Z_required = Phi^-1(1 - (0.05/m)/2)
```

and the minimum detectable excess is `Z_required · sqrt(k·p·(1-p))`, expressed relative
to `k·p`. Because rank 1 is assumed, the resulting figure is a **lower bound** on the
true detection threshold: a gene set not ranked first requires a larger effect.

Values reported in the main text use `k = 100`.

## S3. Per-gene power

Power for a two-sample comparison was computed from the non-central t distribution:

```
power = 1 - F_t(crit; df, ncp) + F_t(-crit; df, ncp)
ncp   = d * sqrt(n1*n2/(n1+n2))
df    = n1 + n2 - 2
```

at α = 0.05 (two-sided) and effect sizes d = 0.5, 0.8, 1.0. These per-gene figures are
reported for context only; because gene-set tests aggregate many genes, the gene-set
detection threshold is better represented by the MDE in §S2.

## S4. Status of the calibration probes

The four cell-type marker sets are excluded from the Benjamini–Hochberg family applied
to candidate gene sets, on the grounds that they are calibration probes rather than
hypotheses under test (main text §2.4). For transparency we report that, had they been
included, none would have reached FDR < 0.05 in any cohort (minimum FDR 0.071 for
microglial markers in GSE213488; 0.217 in GSE256068). Their intended interpretation is
by rank among all tested sets and by observed-versus-expected counts (Table S9).
"""
open(os.path.join(SUP, "Supplementary_Methods.md"), "w", encoding="utf-8").write(sm)

print("[补充材料]")
for f in sorted(os.listdir(SUP)):
    if f.endswith((".csv", ".md")):
        print("   %-52s %7.1f KB" % (f, os.path.getsize(os.path.join(SUP, f)) / 1024))
