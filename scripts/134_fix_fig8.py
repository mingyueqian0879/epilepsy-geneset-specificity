#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 134: 重绘 Fig 8（图例零遮挡）并写回投稿包，同时自检。

做法：直接执行已修好的 scripts/40_fig8_suppl.py，但
  * FIG 重定向到投稿包 04_Figures（这是唯一写回投稿目录的动作）
  * SUP 重定向到临时目录，避免覆盖已交付的补充表；
    同时把临时目录里生成的 S1–S5 与已交付版本逐字节比对，
    作为「输入路径修复是否正确」的独立验证。

自检：任何散点标记的墨迹若触及任一「图例文字」包围盒即报错中止。
"""
import os
import sys
import shutil
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.figure as mfig
from PIL import Image

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = "."
SRC = os.path.join(BASE, "scripts", "40_fig8_suppl.py")
DST = "figures"
SUBDELIV = "results/supplementary"
TMP = os.path.join(BASE, "out", "_fig8_regen")

if os.path.isdir(TMP):
    shutil.rmtree(TMP)
os.makedirs(TMP, exist_ok=True)

CHECK = {"legend_hits": [], "title_clash": [], "fits": []}


def audit(fig, fname):
    if "Fig8" not in os.path.basename(str(fname)):
        return
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for ai, ax in enumerate(fig.axes):
        leg = ax.get_legend()
        if leg is None:
            continue
        boxes = [t.get_window_extent(r) for t in leg.get_texts()]
        pts = []
        for coll in ax.collections:
            offs = coll.get_offsets()
            if offs is None or len(offs) == 0:
                continue
            sizes = coll.get_sizes()
            rad = (np.sqrt(sizes[0]) / 2.0) * (fig.dpi / 72.0) if len(sizes) else 3.0
            pts.append((np.asarray(offs)[:, :2], rad))
        if not pts:
            continue
        xy = np.vstack([p[0] for p in pts])
        rad = max(p[1] for p in pts)
        disp = ax.transData.transform(xy)
        for k, tb in enumerate(boxes):
            dx = np.maximum(np.maximum(tb.x0 - disp[:, 0], 0), disp[:, 0] - tb.x1)
            dy = np.maximum(np.maximum(tb.y0 - disp[:, 1], 0), disp[:, 1] - tb.y1)
            d = np.hypot(dx, dy)
            n = int((d <= rad).sum())
            if n:
                CHECK["legend_hits"].append(
                    "axes%d legend line %d (%r): %d marker(s) within radius %.1f px "
                    "(closest gap %.1f px)"
                    % (ai, k, leg.get_texts()[k].get_text()[:30], n, rad, d.min()))
        fb = leg.get_window_extent(r)
        inv = ax.transAxes.inverted()
        (x0, y0), (x1, y1) = inv.transform([[fb.x0, fb.y0], [fb.x1, fb.y1]])
        if x0 < -0.02 or x1 > 1.0:
            CHECK["fits"].append("axes%d legend x=[%.2f,%.2f] outside axes" % (ai, x0, x1))
        tb = ax.title.get_window_extent(r)
        if not (fb.y1 < tb.y0 or tb.y1 < fb.y0 or fb.x1 < tb.x0 or tb.x1 < fb.x0):
            CHECK["title_clash"].append("axes%d legend overlaps panel title" % ai)


_orig = mfig.Figure.savefig


def patched(self, fname, *a, **kw):
    audit(self, str(fname))
    return _orig(self, fname, *a, **kw)


mfig.Figure.savefig = patched

src = open(SRC, encoding="utf-8").read()
assert 'label_permutation_200perm.csv' in src, "input-path fix missing in 40_fig8_suppl.py"
src = src.replace('FIG = os.path.join(BASE, "figures")', 'FIG = r"%s"' % DST)
src = src.replace('SUP = os.path.join(BASE, "supplementary")', 'SUP = r"%s"' % TMP)
src = src.replace('os.path.join(BASE, "figures")', 'r"%s"' % DST)

print("### regenerating Fig 8 into the submission package")
exec(compile(src, SRC, "exec"), {"__name__": "__main__", "__file__": SRC})

print("\n[self-check]")
ok = True
for cat in ("legend_hits", "title_clash", "fits"):
    if CHECK[cat]:
        ok = False
        for m in CHECK[cat]:
            print("   FAIL [%s] %s" % (cat, m))
    else:
        print("   PASS  no %s" % cat)
if not ok:
    raise SystemExit("FAIL: Fig 8 layout self-check did not pass")

# TIFF 副本
png = os.path.join(DST, "Fig8_reproducibility.png")
im = Image.open(png)
im.convert("RGB").save(os.path.join(DST, "Fig8_reproducibility.tif"), "TIFF",
                       compression="tiff_lzw", dpi=(300, 300))
print("   wrote %s  %dx%d  (%.0f KB)" % (png, im.width, im.height,
                                         os.path.getsize(png) / 1024))

print("\n[consistency] regenerated supplementary tables vs delivered")
for f in sorted(os.listdir(TMP)):
    a = os.path.join(TMP, f)
    b = os.path.join(SUBDELIV, f)
    if not os.path.exists(b):
        print("   %-38s (not in package)" % f)
        continue
    same = open(a, "rb").read() == open(b, "rb").read()
    print("   %-38s %s" % (f, "IDENTICAL" if same else "*** DIFFERS ***"))
