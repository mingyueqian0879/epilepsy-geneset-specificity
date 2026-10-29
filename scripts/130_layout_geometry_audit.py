#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Step 130: 在生成脚本运行时做精确版面几何审计。

做法：劫持 Figure.savefig，在每次落盘前对当前 figure 做：
  (1) 文字包围盒两两重叠（标题/轴标签/刻度/ax.text/图例文字）
  (2) 图例矩形 vs 数据图元（散点集合、柱/矩形、折线）
  (3) 文字越出画布
判定全部在显示坐标系（像素）里做，与 dpi 无关。

注意：这些脚本的输出目录都不是当前投稿包（20260923），运行不会改动投稿文件。
"""
import os
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.figure as mfig
from matplotlib.patches import Rectangle

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

BASE = "."
SCRIPTS = ["37_figures.py", "38_figures2.py", "39_fig5_redo.py",
           "40_fig8_suppl.py", "57_fig7_v2.py", "47_graphabs_v2.py"]

REPORT = {}


def rect_of(artist, renderer):
    bb = artist.get_window_extent(renderer)
    return (bb.x0, bb.y0, bb.x1, bb.y1)


def intersects(a, b):
    return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])


def overlap_area(a, b):
    ix = max(0.0, min(a[2], b[2]) - max(a[0], b[0]))
    iy = max(0.0, min(a[3], b[3]) - max(a[1], b[1]))
    return ix * iy


def text_items(fig, ax, renderer):
    items = []

    def add(tag, t):
        try:
            if t is None:
                return
            s = t.get_text()
            if not s or not s.strip():
                return
            if not t.get_visible():
                return
            items.append((tag, s.strip()[:38], rect_of(t, renderer)))
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
            add("legend", t)
    return items


def data_artists(fig, ax, renderer):
    """返回 [(描述, 矩形列表)]，用于与图例做相交判断。"""
    out = []
    for coll in ax.collections:
        try:
            offs = coll.get_offsets()
            if offs is None or len(offs) == 0:
                continue
            pts = ax.transData.transform(np.asarray(offs)[:, :2])
            sizes = coll.get_sizes()
            rad = (np.sqrt(sizes[0]) / 2.0) * (fig.dpi / 72.0) if len(sizes) else 3.0
            out.append((type(coll).__name__ + "(scatter)", pts, rad))
        except Exception:
            pass
    for p in ax.patches:
        try:
            out.append((type(p).__name__ + "(patch)", [rect_of(p, renderer)], None))
        except Exception:
            pass
    for ln in ax.lines:
        try:
            xy = ln.get_xydata()
            if xy is None or len(xy) == 0:
                continue
            out.append(("Line2D", ax.transData.transform(np.asarray(xy)[:, :2]), None))
        except Exception:
            pass
    return out


def legend_conflicts(fig, ax, renderer):
    leg = ax.get_legend()
    if leg is None:
        return []
    lb = rect_of(leg, renderer)
    hits = []
    for desc, geo, rad in data_artists(fig, ax, renderer):
        if rad is not None:
            # scatter: distance from each marker centre to legend rect
            n = 0
            for (x, y) in geo:
                dx = max(lb[0] - x, 0.0, x - lb[2])
                dy = max(lb[1] - y, 0.0, y - lb[3])
                if np.hypot(dx, dy) <= rad + 0.5:
                    n += 1
            if n:
                hits.append("%s: %d marker(s) within legend rect" % (desc, n))
        elif desc.startswith("Line2D"):
            n = int(((geo[:, 0] >= lb[0]) & (geo[:, 0] <= lb[2])
                     & (geo[:, 1] >= lb[1]) & (geo[:, 1] <= lb[3])).sum())
            if n:
                hits.append("%s: %d vertex/vertices inside legend rect" % (desc, n))
        else:
            for r in geo:
                a = overlap_area(lb, r)
                if a > 0:
                    frac = a / max(1e-9, (r[2] - r[0]) * (r[3] - r[1]))
                    if frac > 0.10:
                        hits.append("%s overlaps legend by %.0f%% of its box"
                                    % (desc, frac * 100))
                        break
    return hits


def audit(fig, fname):
    key = os.path.basename(fname)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.canvas.get_width_height()

    txt_overlaps = []
    outside = []
    leg_hits = []
    for ax in fig.axes:
        items = text_items(fig, ax, r)
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                t1, s1, b1 = items[i]
                t2, s2, b2 = items[j]
                if t1 == "legend" and t2 == "legend":
                    continue
                if not intersects(b1, b2):
                    continue
                a = overlap_area(b1, b2)
                small = min((b1[2] - b1[0]) * (b1[3] - b1[1]),
                            (b2[2] - b2[0]) * (b2[3] - b2[1]))
                if small > 0 and a / small > 0.30:
                    txt_overlaps.append("%s[%s] <-> %s[%s] (%.0f%%)"
                                        % (t1, s1, t2, s2, 100 * a / small))
        for tag, s, b in items:
            # 注意：rcParams 里 savefig.bbox='tight'，落盘时会按所有 artist 的并集
            # 扩张画布，因此 bbox 落在当前 canvas 之外**不等于**被裁切。
            # 这里只在越界幅度极大（> 画布尺寸的 60%）时提示，作为信息项。
            slack_x = 0.6 * W
            slack_y = 0.6 * H
            if (b[0] < -slack_x or b[1] < -slack_y
                    or b[2] > W + slack_x or b[3] > H + slack_y):
                outside.append("%s[%s] bbox=(%.0f,%.0f,%.0f,%.0f) canvas=%dx%d"
                               % (tag, s, b[0], b[1], b[2], b[3], W, H))
        leg_hits.extend(legend_conflicts(fig, ax, r))

    REPORT.setdefault(key, {"text_overlaps": [], "outside": [], "legend": []})
    REPORT[key]["text_overlaps"].extend(txt_overlaps)
    REPORT[key]["outside"].extend(outside)
    REPORT[key]["legend"].extend(leg_hits)


_orig_savefig = mfig.Figure.savefig


def patched_savefig(self, fname, *a, **kw):
    try:
        audit(self, str(fname))
    except Exception as e:
        REPORT.setdefault(os.path.basename(str(fname)), {}).setdefault(
            "audit_errors", []).append(repr(e))
    return _orig_savefig(self, fname, *a, **kw)


mfig.Figure.savefig = patched_savefig

for name in SCRIPTS:
    path = os.path.join(BASE, "scripts", name)
    print("### running %s" % name)
    src = open(path, encoding="utf-8").read()
    try:
        exec(compile(src, path, "exec"), {"__name__": "__main__", "__file__": path})
    except Exception as e:
        print("    SCRIPT FAILED: %r" % e)
    print()

print("#" * 78)
print("LAYOUT AUDIT REPORT")
print("#" * 78)
for k in sorted(REPORT):
    v = REPORT[k]
    probs = (len(v.get("text_overlaps", [])) + len(v.get("outside", []))
             + len(v.get("legend", [])) + len(v.get("audit_errors", [])))
    print("\n%s   %s" % (k, "  <== CHECK" if probs else "  OK"))
    for cat in ("text_overlaps", "legend", "outside", "audit_errors"):
        for msg in v.get(cat, []):
            print("   [%s] %s" % (cat, msg))
