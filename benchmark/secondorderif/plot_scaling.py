#!/usr/bin/env python3
# ============================================================================
# plot_scaling.py — 二阶 Weyl-GTEMPO error-scaling 图（Toulouse / Anderson）
#
# 读取各验证文件夹 results/*.csv，生成 log-log 误差-步长图：
#   secondorderif/toulouse/results/scaling_tau.png / scaling_realtime.png
#   secondorderif/singlemodebath/results/scaling_tau.png / scaling_realtime.png
# ============================================================================
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.abspath(__file__))

CASES = [
    ("toulouse", "Toulouse (U=0, 1 band)",
     [("toulouse_tau.csv", "scaling_tau.png", "Imaginary time ⟨n⟩ (β=1)"),
      ("toulouse_realtime.csv", "scaling_realtime.png", "Real time ⟨n(t=0.4)⟩"),
      ("toulouse_mixed.csv", "scaling_mixed.png", "Mixed contour ⟨n(t=0.4)⟩ (β=1)")]),
    ("singlemodebath", "Anderson (U=0.8, 2 bands)",
     [("anderson_tau.csv", "scaling_tau.png", "Imaginary time ⟨n₁+n₂⟩ (β=1)"),
      ("anderson_realtime.csv", "scaling_realtime.png", "Real time ⟨n₁+n₂(t=0.4)⟩"),
      ("anderson_mixed.csv", "scaling_mixed.png", "Mixed contour ⟨n₁+n₂(t=0.4)⟩ (β=1)")]),
]

COLORS = {1: "tab:blue", 2: "tab:red"}


def read_csv(path):
    pts = defaultdict(list)   # order -> (h, err)
    with open(path) as f:
        for row in csv.DictReader(f):
            pts[int(row["order"])].append((float(row["h"]), float(row["abs_err"])))
    return {o: sorted(v) for o, v in pts.items()}


def slope(hs, errs):
    import math
    xs = [math.log(h) for h in hs]
    ys = [math.log(e) for e in errs]
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


for folder, title, jobs in CASES:
    res = os.path.join(REPO, folder, "results")
    for csvname, pngname, ylabel in jobs:
        path = os.path.join(res, csvname)
        if not os.path.exists(path):
            print("skip (missing):", path)
            continue
        data = read_csv(path)
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
        for order in sorted(data):
            hs = [p[0] for p in data[order]]
            errs = [p[1] for p in data[order]]
            mk = "o-" if order == 2 else "s--"
            ax.loglog(hs, errs, mk, color=COLORS[order],
                      label=f"order {order}  (slope {slope(hs, errs):.3f})")
        # O(h²) 参考线
        ref = data[2]
        h0, e0 = ref[0]
        hs_ref = [h0, h0 / 4]
        ax.loglog(hs_ref, [e0, e0 / 16], "k:", lw=1, label=r"$\propto h^2$")
        ax.set_xlabel("step  h")
        ax.set_ylabel("|error|")
        ax.set_title(f"{ylabel}\n{title} — Weyl 2nd order vs ED", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        out = os.path.join(res, pngname)
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print("written:", out)
