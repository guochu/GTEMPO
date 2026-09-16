# Plot all IDMRG1 kvals (per-site residuals of every sweep) into debug/mult/kvals.png
# Run: /tmp/triqs_env/bin/python debug/mult/plot_kvals.py
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))

cases = sorted(glob.glob(os.path.join(here, "kvals_*_D*.csv")))
ncase = len(cases)
ntrial = 3
fig, axes = plt.subplots(ntrial, ncase, figsize=(6.5 * ncase, 3.6 * ntrial),
                         squeeze=False, sharex=False)

for c, csv in enumerate(cases):
    data = {}          # (trial, sweep) -> (seg, [residuals in sweep order])
    for line in open(csv):
        parts = line.strip().split(",")
        if parts[0] == "trial":
            continue
        trial, sweep, seg, pos = int(parts[0]), int(parts[1]), parts[2], int(parts[3])
        res = float(parts[5])
        data.setdefault((trial, sweep), [[], []])
        (data[(trial, sweep)][0] if seg == "left" else data[(trial, sweep)][1])
        data[(trial, sweep)][0 if seg == "left" else 1].append((pos, res))

    tag = os.path.basename(csv)[len("kvals_"):-len(".csv")]
    for (trial, sweep), (left, right) in sorted(data.items()):
        left.sort()
        right.sort(key=lambda t: -t[0])     # right sweep runs sites L..2
        seq = [r for _, r in left] + [r for _, r in right]
        ax = axes[trial - 1][c]
        xs = np.arange(1, len(seq) + 1)
        ax.semilogy(xs, seq, "-", lw=1, color=plt.cm.viridis(sweep / 12), label=f"sweep {sweep}")
        # mark the left->right junction
        if left:
            ax.axvline(len(left) + 0.5, color="gray", ls=":", lw=0.8)
    ax = axes[0][c]
    ax.set_title(tag.replace("_", " "), fontsize=11)
    ax.set_ylabel("kvals (residual norm)")
    axes[-1][c].set_xlabel("site index within sweep (left 1..L-1 | right L..2)")
    for r in range(ntrial):
        axes[r][c].set_ylim(bottom=max(1e-16, axes[r][c].get_ylim()[0]))

handles, labels = axes[0][0].get_legend_handles_labels()
fig.legend(handles, labels, fontsize=8, loc="upper center", ncol=min(12, len(labels)),
           bbox_to_anchor=(0.5, 0.995))
fig.suptitle("IDMRG1 mult: per-site kvals of every sweep "
             "(dotted line = left/right junction)", fontsize=12, y=1.0)
fig.tight_layout(rect=(0, 0, 1, 0.90))
out = os.path.join(here, "kvals.png")
fig.savefig(out, dpi=130)
print("written", out)
