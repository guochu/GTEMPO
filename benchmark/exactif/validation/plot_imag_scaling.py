# plot_imag_scaling.py —— 虚轴 Trotter error scaling 干净图
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 0.8,
    "figure.dpi": 160,
})

rows = list(csv.DictReader(open("imag_scaling.csv")))
data = {}
for r in rows:
    data.setdefault(int(r["order"]), []).append((int(r["N"]), float(r["nerr"])))
for k in data:
    data[k].sort()

C = {1: "#1f77b4", 2: "#d62728", 4: "#2ca02c"}
MK = {1: "o", 2: "s", 4: "D"}
LB = {1: "order 1", 2: "order 2", 4: "order 4"}

fig, ax = plt.subplots(figsize=(7.2, 5.0))

# 阶梯参考线 p=2 / p=4
ns = [1, 2, 4, 8, 16, 32, 64]
e2 = [data[1][0][1] * (n / data[1][0][0]) ** -2.0 for n in ns]
e4 = [data[4][0][1] * (n / data[4][0][0]) ** -4.0 for n in ns]
ax.loglog(ns, e2, ":", color="gray", lw=1.0, label=r"$\propto \Delta\tau^2$")
ax.loglog(ns, e4, "--", color="gray", lw=1.0, label=r"$\propto \Delta\tau^4$")

for order in (1, 2, 4):
    xs = [d[0] for d in data[order]]
    ys = [d[1] for d in data[order]]
    ms = 9 if order != 4 else 7
    ax.loglog(xs, ys, marker=MK[order], color=C[order], lw=1.4,
              markersize=ms, markerfacecolor="white", markeredgewidth=1.3,
              label=LB[order])

# 数值斜率标注（最后三点拟合）
import math
def slope(k):
    pts = data[k][-3:]
    x0 = math.log(pts[0][0]); x1 = math.log(pts[-1][0])
    y0 = math.log(pts[0][1]); y1 = math.log(pts[-1][1])
    return (y1 - y0) / (x1 - x0)
ax.annotate(f"convergence  p = {abs(slope(1)):.2f}", xy=(22, data[1][-2][1]), color=C[1],
            fontsize=10, xytext=(26, data[1][-2][1] * 1.35))
ax.annotate(f"convergence  p = {abs(slope(4)):.2f}", xy=(22, data[4][-2][1]), color=C[4],
            fontsize=10, xytext=(26, data[4][-2][1] * 0.6))

ax.set_xlabel(r"$N$  ($\Delta\tau=\beta/N$)", fontsize=11)
ax.set_ylabel(r"impurity density error  $|n(\epsilon_d)-n_\ast|$", fontsize=11)
ax.set_title("Imaginary-axis strict-IF: Trotter error scaling"
             "  ($\\beta$=1, $\\epsilon_d$=0.7, $W_B$=1, $A_B$=0.5)", fontsize=11.5)
ax.set_xticks(ns)
ax.set_xticklabels([str(n) for n in ns])
ax.tick_params(which="both", direction="in", top=True, right=True)
ax.grid(True, which="major", ls=":", lw=0.5, alpha=0.6)
ax.legend(loc="lower left", fontsize=10, frameon=True, framealpha=0.9)
ax.annotate("orders 1 and 2 coincide on the trace loop\n"
            "(trace cyclicity cancels odd-order terms)",
            xy=(2.2, 3.4e-4), fontsize=8.5, color="#555555")
ax.annotate("GMPS single-MPS pipeline reproduces the\n"
            "discrete chain to machine precision (≤1.3e-13)\n"
            "→ error shown = schedule Trotter error",
            xy=(1.15, 5.2e-3), fontsize=8.5, color="#555555",
            xytext=(1.15, 2.2e-5),
            arrowprops=dict(arrowstyle="->", color="#555555", lw=0.8))

plt.tight_layout()
plt.savefig("imag_scaling.png", bbox_inches="tight")
print("saved imag_scaling.png")
