# Run with a python that has numpy + matplotlib (e.g. the TRIQS env):
#   python docs/tutorials/bethedmft/imag/compare_plot.py
#
# Compares the Bethe-DMFT iteration sequences produced by the two impurity
# solvers — GTEMPO (imaginary-time Grassmann tensor networks) and TRIQS
# cthyb (CT-QMC) — which are fed with the identical per-iteration
# hybridisations Δ(iω) stored in data/gtempo_b*.json.
#
# For each U the script draws one figure with six panels showing the G(τ)
# comparison at DMFT iterations 1, 3, 5, 6, 7, 8, and reports the
# per-iteration deviations (also written to data/consistency.txt).

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, "data")

sets = [(5.0, 1.0), (5.0, 5.0), (10.0, 1.0), (10.0, 5.0)]
iters_shown = [1, 2, 3, 4, 5, 6]          # iterations shown in the panels
report = []


def load(solver, tag):
    path = os.path.join(datadir, f"{solver}_{tag}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


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
        g_g = np.array(itg["G_tau"])
        tau_q = np.array(itq["tau"])
        g_q = np.interp(tau_g, tau_q, np.array(itq["G_tau"]))
        err_tau = np.linalg.norm(g_g - g_q) / np.linalg.norm(g_g)
        # G(iω): match the two solvers' frequency grids (positive ωₙ only)
        ws_g = np.array(itg["ws"])
        ws_q = np.array(itq["iw"])[:, 1]
        lookup = {round(w, 9): idx for idx, w in enumerate(ws_g)}
        pairs = [(lookup[round(w, 9)], j) for j, w in enumerate(ws_q)
                 if round(w, 9) in lookup]
        ig = [p[0] for p in pairs]
        iq = [p[1] for p in pairs]
        giw_all_g = np.array(itg["G_iw"])[:, 0] + 1j * np.array(itg["G_iw"])[:, 1]
        giw_all_q = np.array(itq["G_iw"])[:, 0] + 1j * np.array(itq["G_iw"])[:, 1]
        giw_g, giw_q = giw_all_g[ig], giw_all_q[iq]
        err_iw = np.linalg.norm(giw_g - giw_q) / np.linalg.norm(giw_q)
        lines.append((k, err_tau, err_iw))
        print(f"  iter {k:2d}: rel.dev G(τ) = {err_tau:.3e}   "
              f"rel.dev G(iω) = {err_iw:.3e}   ({len(pairs)} freqs)")
    report.append((beta, U, lines))

    # ---------- G(τ) panels at selected iterations ----------
    avail = [k for k in iters_shown if k in by_iter_g]
    ncol = 3
    nrow = (len(avail) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(13, 3.4 * nrow),
                             sharex=True, sharey=True)
    for ax, k in zip(np.atleast_1d(axes).ravel(), avail):
        itg, itq = by_iter_g[k], by_iter_q[k]
        tau_g = np.array(itg["tau"])
        tau_q = np.array(itq["tau"])
        g_q = np.interp(tau_g, tau_q, np.array(itq["G_tau"]))
        ax.plot(tau_g, itg["G_tau"], "-", lw=2, label="GTEMPO")
        ax.plot(tau_g, g_q, "o", ms=2.5, mfc="none", label="cthyb (CT-QMC)")
        err = np.linalg.norm(np.array(itg["G_tau"]) - g_q) \
            / np.linalg.norm(np.array(itg["G_tau"]))
        ax.set_title(f"iteration {k}   (rel. dev. {err:.1e})", fontsize=10)
        ax.set_xlabel("τ"); ax.set_ylabel("G(τ)")
    axes.flat[0].legend(fontsize=8)
    fig.suptitle(f"β = {beta}, U = {U}: G(τ) per DMFT iteration, "
                 "identical Δ(iω) sequence", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(here, f"compare_gtau_{tag}.png"), dpi=130)
    plt.close(fig)
    print(f"written compare_gtau_{tag}.png")

with open(os.path.join(datadir, "consistency.txt"), "w") as f:
    for (beta, U, lines) in report:
        f.write(f"beta = {beta}, U = {U}\n")
        f.write("iter   rel.dev G(tau)   rel.dev G(iw)\n")
        for (k, e1, e2) in lines:
            f.write(f"{k:4d}   {e1:.6e}   {e2:.6e}\n")
        f.write("\n")
print("\nwritten consistency.txt")
