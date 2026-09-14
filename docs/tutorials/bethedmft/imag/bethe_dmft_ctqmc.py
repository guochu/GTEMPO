# Run with the TRIQS environment python:
#   python docs/tutorials/bethedmft/imag/bethe_dmft_ctqmc.py [set]
#
# CT-QMC (TRIQS cthyb) solution of the same Bethe-DMFT iteration sequence
# as bethe_dmft_gtempo.jl: the per-iteration Matsubara hybridisations are
# read from the pole parameters stored in data/gtempo_b*.json, so both
# solvers are fed with *identical* Δ(iω) sequences and start from the
# *identical* initial state.  Results are written to data/ctqmc_b*.json.
#
# Optional env vars:
#   DMFT_SET=k        run only parameter set k (1..4), same indexing as julia
#   DMFT_NCYCLES=n    QMC cycles per iteration (default 200000)

import json
import os
import time
import numpy as np
from triqs.gf import Gf, MeshImFreq, MeshImTime, iOmega_n, inverse
from triqs.operators import n
from triqs_cthyb import Solver

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, "data")
set_idx = int(os.environ.get("DMFT_SET", "0"))
n_cycles = int(os.environ.get("DMFT_NCYCLES", "200000"))
seed = 777

sets = [(5.0, 1.0), (5.0, 5.0), (10.0, 1.0), (10.0, 5.0)]
run_sets = sets if set_idx == 0 else [sets[set_idx - 1]]

for (beta, U) in run_sets:
    tag = f"b{int(beta)}_U{int(U)}"
    src = os.path.join(datadir, f"gtempo_{tag}.json")
    if not os.path.exists(src):
        print(f"skip β = {beta}, U = {U}: {src} not found", flush=True)
        continue
    with open(src) as f:
        gtempo = json.load(f)
    if "ws" not in gtempo["iterations"][0]:
        print(f"skip β = {beta}, U = {U}: stale file without Δ(iω) values "
              "(rerun bethe_dmft_gtempo.jl)", flush=True)
        continue
    nw, ntau = gtempo["nw"], gtempo["ntau"]
    mu_imp = gtempo["mu_imp"]          # ϵ_d = −U/2 (GTEMPO level convention)
    mu_ctqmc = -mu_imp                 # cthyb: G₀⁻¹ = iω + mu − Δ  ⇒  mu = U/2

    print(f"=== Bethe DMFT (cthyb): β = {beta}, U = {U}, n_iw = {nw}, "
          f"n_cycles = {n_cycles} ===", flush=True)

    iterations = []
    for itinfo in gtempo["iterations"]:
        it = itinfo["iter"]
        # Δ(iω) values from the GTEMPO export, looked up by ωₙ = Im(iωₙ)
        ws_all = np.array(itinfo["ws"])
        delta_all = np.array(itinfo["Delta_iw"])[:, 0] + \
            1j * np.array(itinfo["Delta_iw"])[:, 1]
        dmap = {round(w, 9): d for w, d in zip(ws_all, delta_all)}
        # the cthyb tail-fitter needs a minimum number of Matsubara points;
        # only the first nw positive frequencies are compared with GTEMPO
        n_iw_solver = max(nw, 40)
        n_tau_solver = max(ntau + 1, 6 * n_iw_solver + 2)
        S = Solver(beta=beta, n_iw=n_iw_solver, n_tau=n_tau_solver,
                   gf_struct=[("up", 1), ("down", 1)])
        iw = np.array([complex(w) for w in S.G0_iw["up"].mesh])   # iωₙ values
        Delta = S.G0_iw["up"].copy()
        Delta.zero()
        for k, w in enumerate(iw):
            Delta.data[k, 0, 0] = dmap[round(w.imag, 9)]
        for spin in ("up", "down"):
            S.G0_iw[spin] << inverse(iOmega_n + mu_ctqmc - Delta)

        t0 = time.time()
        S.solve(h_int=U * n("up", 0) * n("down", 0),
                n_cycles=n_cycles, length_cycle=max(50, int(beta * 10)),
                n_warmup_cycles=2000, random_seed=seed)
        wall = time.time() - t0

        G_iw = 0.5 * (S.G_iw["up"] + S.G_iw["down"])
        G_tau = 0.5 * (S.G_tau["up"] + S.G_tau["down"])
        giw_all = np.array(G_iw.data[:, 0, 0])
        gtau = np.array(G_tau.data[:, 0, 0])
        tau = np.array(list(G_tau.mesh))
        # positive Matsubara frequencies n = 0..nw-1 (GTEMPO grid)
        iw_pos = np.array([complex(w) for w in G_iw.mesh])
        order = sorted(range(len(iw_pos)), key=lambda k: iw_pos[k].imag)
        pos = [k for k in order if iw_pos[k].imag > 0][:nw]
        iterations.append({
            "iter": it,
            "iw": [[float(w.real), float(w.imag)] for w in (iw_pos[k] for k in pos)],
            "G_iw": [[float(giw_all[k].real), float(giw_all[k].imag)] for k in pos],
            "G_tau": [float(g) for g in gtau],
            "tau": [float(t) for t in tau],
        })
        print(f"  iter {it:2d} done ({wall:.1f} s, "
              f"sign = {S.average_sign:.3f})", flush=True)

    out = {"solver": "cthyb", "beta": beta, "U": U, "nw": nw, "ntau": ntau,
           "n_cycles": n_cycles, "iterations": iterations}
    with open(os.path.join(datadir, f"ctqmc_{tag}.json"), "w") as f:
        json.dump(out, f)
    print(f"written ctqmc_{tag}.json", flush=True)

print("all done")
