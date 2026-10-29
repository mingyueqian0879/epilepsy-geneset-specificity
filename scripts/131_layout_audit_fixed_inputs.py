#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 131: 修正输入路径后的版面几何审计（不触碰任何交付目录）。

对每个生成脚本做两处替换后再执行：
  (1) 输入改用 label_permutation_200perm.csv
      —— label_permutation.csv 现已被 100 次置换的结果覆盖（与
      label_permutation_100perm.csv 逐字节相同），而图 5 与表 S5 需要 200 次
      置换的列（null_median / null_p95 / emp_P），脚本因此全部报错。
  (2) 所有输出目录重定向到 out/_layout_repro/，交付目录只读。

savefig 被劫持：落盘前审计「文字压文字」「图例压数据」，并给出被压点的数据坐标。
"""
import os
import sys
import shutil
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.figure as mfig

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = "."
TMP = os.path.join(BASE, "out", "_layout_repro")
SUB = "figures"
os.makedirs(TMP, exist_ok=True)

JOBS = [
    ("37_figures.py", "Fig1_field_survey.png"),
    # 38 先跑：它也会写一份 Fig5/Fig7（旧版），随后被 39/57 覆盖为交付版
    ("38_figures2.py", "Fig6_cell_composition.png"),
    ("39_fig5_redo.py", "Fig5_label_permutation.png"),
    ("40_fig8_suppl.py", "Fig8_reproducibility.png"),
    ("57_fig7_v2.py", "Fig7_reporting_standard.png"),
    ("47_graphabs_v2.py", "graphical_abstract.png"),
    ("66_figs_update.py", "Fig3_Z_distribution.png"),
    ("66_figs_update.py", "Fig4_hot_vs_positive.png"),
    ("66_figs_update.py", "Fig9_spikein_validation.png"),
]

REPORT = {}


def rect_of(a, r):
    b = a.get_window_extent(r)
    return (b.x0, b.y0, b.x1, b.y1)


def inter_area(a, b):
    ix = max(0.0, min(a[2], b[2]) - max(a[0], b[0]))
    iy = max(0.0, min(a[3], b[3]) - max(a[1], b[1]))
    return ix * iy


def audit(fig, fname):
    key = os.path.basename(str(fname))
    rec = REPORT.setdefault(key, {"text": [], "legend": [], "errors": []})
    try:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        for ai, ax in enumerate(fig.axes):
            items = []

            def add(tag, t):
                try:
                    if t is None or not t.get_visible():
                        return
                    s = (t.get_text() or "").strip()
                    if s:
                        items.append((tag, s[:34].replace("\n", "|"), rect_of(t, r)))
                except Exception:
                    pass

            add("title", ax.title)
            add("xlabel", ax.xaxis.label)
            add("ylabel", ax.yaxis.label)
            for t in ax.get_xticklabels():
                add("xtick", t)
            for t in ax.get_yticklabels():
                add("ytick", t)
            for t in ax.texts:
                add("axtext", t)
            leg = ax.get_legend()
            if leg is not None:
                for t in leg.get_texts():
                    add("legtext", t)

            for i in range(len(items)):
                for j in range(i + 1, len(items)):
                    t1, s1, b1 = items[i]
                    t2, s2, b2 = items[j]
                    if t1 == "legtext" and t2 == "legtext":
                        continue
                    if not (b1[2] <= b2[0] or b2[2] <= b1[0]
                            or b1[3] <= b2[1] or b2[3] <= b1[1]):
                        a = inter_area(b1, b2)
                        small = min((b1[2] - b1[0]) * (b1[3] - b1[1]),
                                    (b2[2] - b2[0]) * (b2[3] - b2[1]))
                        if small > 0 and a / small > 0.30:
                            rec["text"].append("axes%d %s[%s] <-> %s[%s] %.0f%%"
                                               % (ai, t1, s1, t2, s2, 100 * a / small))

            # legend vs data
            if leg is not None:
                lb = rect_of(leg, r)
                inv = ax.transAxes.inverted()
                (lx0, ly0), (lx1, ly1) = inv.transform([[lb[0], lb[1]], [lb[2], lb[3]]])
                for coll in ax.collections:
                    offs = coll.get_offsets()
                    if offs is None or len(offs) == 0:
                        continue
                    data = np.asarray(offs)[:, :2]
                    disp = ax.transData.transform(data)
                    sizes = coll.get_sizes()
                    rad = (np.sqrt(sizes[0]) / 2.0) * (fig.dpi / 72.0) if len(sizes) else 3.0
                    dx = np.maximum(np.maximum(lb[0] - disp[:, 0], 0), disp[:, 0] - lb[2])
                    dy = np.maximum(np.maximum(lb[1] - disp[:, 1], 0), disp[:, 1] - lb[3])
                    d = np.hypot(dx, dy)
                    idx = np.where(d <= rad + 0.5)[0]
                    for k in idx:
                        rec["legend"].append(
                            "axes%d LEGEND covers a marker: data=(%.2f, %.2f) "
                            "legend x=[%.2f,%.2f] y=[%.2f,%.2f]"
                            % (ai, data[k, 0], data[k, 1], lx0, lx1, ly0, ly1))
                for p in ax.patches:
                    pb = rect_of(p, r)
                    a = inter_area(lb, pb)
                    if a > 0:
                        frac = a / max(1e-9, (pb[2] - pb[0]) * (pb[3] - pb[1]))
                        if frac > 0.10:
                            rec["legend"].append(
                                "axes%d LEGEND covers a %s (%.0f%% of its box)"
                                % (ai, type(p).__name__, 100 * frac))
                for ln in ax.lines:
                    xy = ln.get_xydata()
                    if xy is None or len(xy) == 0:
                        continue
                    d2 = ax.transData.transform(np.asarray(xy)[:, :2])
                    n = int(((d2[:, 0] >= lb[0]) & (d2[:, 0] <= lb[2])
                             & (d2[:, 1] >= lb[1]) & (d2[:, 1] <= lb[3])).sum())
                    if n:
                        rec["legend"].append("axes%d LEGEND covers %d line vertices"
                                             % (ai, n))
    except Exception as e:
        rec["errors"].append(repr(e))


_orig = mfig.Figure.savefig


def patched(self, fname, *a, **kw):
    audit(self, str(fname))
    return _orig(self, fname, *a, **kw)


mfig.Figure.savefig = patched

# 只跑需要的脚本一次，取其中关心的图
done = set()
for script, target in JOBS:
    if script in done:
        continue
    done.add(script)
    path = os.path.join(BASE, "scripts", script)
    src = open(path, encoding="utf-8").read()
    src = src.replace('"label_permutation.csv"', '"label_permutation_200perm.csv"')
    src = src.replace('os.path.join(BASE, "figures")', 'r"%s"' % TMP)
    src = src.replace('os.path.join(BASE, "supplementary")', 'r"%s"' % TMP)
    src = src.replace(r'"04_Figures"', 'r"%s"' % TMP)
    print("### %s" % script)
    try:
        exec(compile(src, path, "exec"), {"__name__": "__main__", "__file__": path})
    except Exception as e:
        print("    FAILED: %r" % e)
    print()

print("#" * 78)
print("GEOMETRY AUDIT (corrected inputs, outputs redirected)")
print("#" * 78)
for script, target in JOBS:
    rec = REPORT.get(target)
    if rec is None:
        print("\n%-34s NOT PRODUCED" % target)
        continue
    n = len(rec["text"]) + len(rec["legend"]) + len(rec["errors"])
    print("\n%-34s %s" % (target, "  <== CHECK" if n else "  OK"))
    for cat in ("legend", "text", "errors"):
        seen = set()
        for m in rec[cat]:
            if m in seen:
                continue
            seen.add(m)
            print("   [%s] %s" % (cat, m))

print("\n[reproduced outputs in] %s" % TMP)
for f in sorted(os.listdir(TMP)):
    if f.endswith(".png"):
        a = os.path.join(TMP, f)
        b = os.path.join(SUB, f)
        print("   %-34s reproduced=%8d  delivered=%s"
              % (f, os.path.getsize(a),
                 os.path.getsize(b) if os.path.exists(b) else "ABSENT"))
