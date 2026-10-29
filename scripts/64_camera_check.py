#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 64: 决定性问题 —— camera 找到的 22 个"显著"基因集是否生物学合理?
  若同样是不合理的(味觉转导/寄生虫感染等) -> 说明本文方法校准更好
  若是合理的通路 -> 说明本文方法过保守, 结论须软化
"""
import pandas as pd
import numpy as np

OUT = "out"
C = pd.read_csv(OUT + r"\method_comparison.csv")

for ds in C.dataset.unique():
    s = C[C.dataset == ds].copy()
    print("=" * 94)
    print(" [%s]  camera 检出的基因集 (FDR<0.05)" % ds)
    print("=" * 94)
    cam = s[s.FDR_camera < 0.05].sort_values("FDR_camera")
    if not len(cam):
        print("  (无)")
    else:
        print("  %-48s %-6s %-9s %s" % ("基因集", "n", "本文FDR", "camera FDR"))
        for _, r in cam.iterrows():
            print("  %-48s %-6d %-9.3f %.4f" % (r.gene_set[:48], r.n,
                  r.FDR_matched if pd.notna(r.FDR_matched) else np.nan, r.FDR_camera))
    print()

    # 生物学合理性评分: 检查是否含明显无关主题
    IRRELEVANT = ["TASTE", "LEISHMANIA", "PRION", "MYOGENESIS", "ADIPOGENESIS",
                  "BILE_ACID", "TERPENOID", "BETA_ALANINE", "GLYCOSPHINGOLIPID",
                  "FATTY_ACID", "STARCH", "SUCROSE", "HEDGEHOG", "DILATED_CARDIOMYOPATHY",
                  "ASTHMA", "MALARIA", "AMOEBIASIS", "CHAGAS", "LEGIONELLOSIS",
                  "NTSE", "PROSTATE_CANCER", "BLADDER_CANCER", "THYROID_CANCER",
                  "ENDOMETRIAL_CANCER", "GLIOMA", "MELANOMA", "COLORECTAL"]
    cam_irr = cam.gene_set.str.upper().str.contains("|".join(IRRELEVANT)).sum()
    print("  camera 显著集中, 明显与癫痫无关的主题 = %d / %d (%.0f%%)" % (
        cam_irr, len(cam), 100 * cam_irr / max(len(cam), 1)))

    # 本文方法排前 10 的集
    print("\n  本文方法 Z 最高的 10 个集 (供对比):")
    top = s.dropna(subset=["Z_matched"]).nlargest(10, "Z_matched")
    for _, r in top.iterrows():
        mark = "  <- 不合理主题" if any(k in r.gene_set.upper() for k in IRRELEVANT) else ""
        print("     %-48s Z=%+.2f  FDR_本文=%.3f  FDR_camera=%.3f%s" % (
            r.gene_set[:48], r.Z_matched, r.FDR_matched,
            r.FDR_camera if pd.notna(r.FDR_camera) else np.nan, mark))
    print()

    # camera 与本文的 TOP 集合重叠
    cam_top = set(cam.gene_set.head(20))
    my_top = set(s.dropna(subset=["Z_matched"]).nlargest(20, "Z_matched").gene_set)
    print("  camera TOP20 与 本文 TOP20 重叠 = %d 个" % len(cam_top & my_top))
    if cam_top & my_top:
        print("     重叠: %s" % ", ".join(list(cam_top & my_top)[:8]))
    print()
