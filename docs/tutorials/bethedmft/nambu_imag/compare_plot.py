# Run with a python that has numpy + matplotlib (e.g. the TRIQS env):
#   python docs/tutorials/bethedmft/nambu_imag/compare_plot.py
#
# Compares the superconducting (Nambu) Bethe-DMFT iteration sequences
# produced by the two impurity solvers — GTEMPO (imaginary-time Grassmann
# tensor networks) and TRIQS cthyb (CT-QMC) — which are fed with the
# identical per-iteration Nambu hybridisations stored in
# data/gtempo_nambu_b*.json.
#
# cthyb is run in the Werner–Millis matrix mode with abstract flavours
# c₀ = d↑, c₁ = d↓†; its raw measurements are translated to the physical
# components via
#     G↑↑ = G₀₀ ,  G↓↓ = G₁₁(β−τ) ,  ⟨T d↑d↓⟩ = −G₀₁ ,  ⟨T d↓d↑⟩ = −G₁₀ .
#
# For each U the script draws one figure with six panels (DMFT iterations
# 1..6); each panel shows the four physical G components, GTEMPO as lines
# and cthyb as symbols (the anomalous components on a magnified twin axis),
# and reports the per-iteration deviations (data/consistency_nambu.txt).

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, "data")

sets = [(1.0, 1.0), (1.0, 5.0)]
iters_shown = [1, 2, 3, 4, 5, 6]
comps = ["uu", "dd", "ud", "du"]
comp_labels = {"uu": r"$G_{uu}$", "dd": r"$G_{dd}$",
               "ud": r"$G_{ud}$", "du": r"$G_{du}$"}
anom = {"ud", "du"}
report = []


def load(solver, tag):
    path = os.path.join(datadir, f"{solver}_nambu_{tag}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def physical_cthyb(itq):
    """abstract cthyb G_tau -> physical components on the cthyb tau mesh"""
    tau = np.array(itq["tau"])
    return {
        "uu": np.array(itq["G_tau_00"]),
        "dd": np.array(itq["G_tau_11"])[::-1],
        "ud": -np.array(itq["G_tau_01"]),
        "du": -np.array(itq["G_tau_10"])[::-1],
    }


for (beta, U) in sets:
    tag = f"b{int(beta)}_U{int(U)}"
    gt, qmc = load("gtempo", tag), load("ctqmc", tag)
    if gt is None or qmc is None:
        print(f"skip β = {beta}, U = {U}: missing solver data")
        continue
    print(f"\n=== β = {beta}, U = {U} ===")

    by_iter_g = {it["iter"]: it for it in gt["iterations"]}
    by_iter_q = {it["iter"]: it for it in qmc["iterations"]}

    # ---------- per-iteration consistency ----------
    lines = []
    for k in sorted(by_iter_g):
        itg, itq = by_iter_g[k], by_iter_q[k]
        tau_g = np.array(itg["tau"])
        phys_q = physical_cthyb(itq)
        errs = {}
        for c in comps:
            g_g = np.array(itg[f"G_tau_{c}"])
            g_q = np.interp(tau_g, np.array(itq["tau"]), phys_q[c])
            errs[c] = np.linalg.norm(g_g - g_q) / np.linalg.norm(g_g)
        lines.append((k, errs))
        print("  iter %2d: " % k
              + "  ".join(f"{c}: {errs[c]:.3e}" for c in comps))
    report.append((beta, U, lines))

    # ---------- G(τ) panels at selected iterations ----------
    avail = [k for k in iters_shown if k in by_iter_g]
    ncol = 3
    nrow = (len(avail) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(13, 4.2 * nrow),
                             sharex=True)
    legend_handles = None
    for ax, k in zip(np.atleast_1d(axes).ravel(), avail):
        itg, itq = by_iter_g[k], by_iter_q[k]
        tau_g = np.array(itg["tau"])
        tau_q = np.array(itq["tau"])
        phys_q = physical_cthyb(itq)
        axt = ax.twinx()
        errs = {}
        hs, ls = [], []
        for c in comps:
            g_g = np.array(itg[f"G_tau_{c}"])
            g_q = np.interp(tau_g, tau_q, phys_q[c])
            errs[c] = np.linalg.norm(g_g - g_q) / np.linalg.norm(g_g)
            color = f"C{comps.index(c)}"
            target = axt if c in anom else ax
            ls_ = "--" if c in anom else "-"
            hs.append(target.plot(tau_g, g_g, ls_, lw=2, color=color,
                                  label=f"GTEMPO {comp_labels[c]}")[0])
            hs.append(target.plot(tau_g, g_q, "o", ms=2.5, mfc="none",
                                  color=color,
                                  label=f"cthyb {comp_labels[c]}")[0])
            ls.append(f"GTEMPO {comp_labels[c]}")
            ls.append(f"cthyb {comp_labels[c]}")
        if legend_handles is None:
            legend_handles = (hs, ls)
        emax = max(errs.values())
        ax.set_title(f"iteration {k}   (max rel. dev. {emax:.1e})",
                     fontsize=10)
        ax.set_xlabel("τ")
        ax.set_ylabel("normal  $G_{uu},G_{dd}$")
        axt.set_ylabel("anomalous  $G_{ud},G_{du}$")
        axt.tick_params(axis="y", colors="gray")
    fig.legend(legend_handles[0], legend_handles[1], fontsize=8,
               loc="upper center", ncol=4, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(f"β = {beta}, U = {U}: Nambu DMFT G(τ) per iteration, "
                 "identical Δ(iω) sequence (anomalous components dashed, "
                 "right axis)", fontsize=12, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(os.path.join(here, f"compare_gtau_{tag}.png"), dpi=130)
    plt.close(fig)
    print(f"written compare_gtau_{tag}.png")

with open(os.path.join(datadir, "consistency_nambu.txt"), "w") as f:
    for (beta, U, lines) in report:
        f.write(f"beta = {beta}, U = {U}\n")
        f.write("iter   " + "   ".join(f"rel.dev G(tau) {c}" for c in comps)
                + "\n")
        for (k, errs) in lines:
            f.write(f"{k:4d}   " + "   ".join(f"{errs[c]:.6e}"
                                              for c in comps) + "\n")
        f.write("\n")
print("\nwritten consistency_nambu.txt")
