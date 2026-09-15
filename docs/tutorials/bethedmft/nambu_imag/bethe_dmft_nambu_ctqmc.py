# Run with the TRIQS environment python:
#   python docs/tutorials/bethedmft/nambu_imag/bethe_dmft_nambu_ctqmc.py
#
# CT-QMC (TRIQS cthyb) solution of the same superconducting (Nambu) Bethe-DMFT
# iteration sequence as bethe_dmft_nambu_gtempo.jl.  The per-iteration
# Matsubara hybridisations (four blocks: uu, dd normal; ud, du anomalous) are
# read from data/gtempo_nambu_b*.json, so both solvers are fed with the
# *identical* Δ(iω) sequence and start from the *identical* initial state.
#
# Nambu formulation (Werner–Millis matrix mode, delta interface): one block
# 'nm' of size 2 whose flavors are the quasiparticles c₀ = d↑ and c₁ = d↓†.
# Inputs:
#   * S.Delta_tau['nm'] : the 2×2 hybridisation matrix on the τ mesh.  With the
#     abstract-flavour action S_hyb = ∫∫ c̄ᵢ(τ) Δᵢⱼ(τ−τ′) cⱼ(τ′) the four GTEMPO
#     blocks enter *directly* (no p-h conjugation!):
#       Δ₀₀ = Δuu (normal ↑),  Δ₁₁ = Δdd (hole-channel kernel of the ↓ leg),
#       Δ₀₁ = Δud (pair creation),  Δ₁₀ = Δdu (pair annihilation).
#   * h_loc0 = ε_d (n₀ − n₁) with ε_d = −U/2 passed to solve() (the hole-like
#     flavor carries the opposite on-site level; n₁ = 1 − n↓),
#   * h_int = U n↑n↓ = U n₀ (1 − n₁).
# The measured abstract G's translate to the physical components as
#   G↑↑ = G₀₀ ,  G↓↓ = −G₁₁(β−τ) ,  ⟨T d↑d↓⟩ = −G₀₁ ,  ⟨T d↓d↑⟩ = −G₁₀ .
# The off-diagonal (anomalous) blocks make cthyb run in its generic Matrix
# mode (no segment picture).
#
# Optional env vars:
#   DMFT_SET=k        run only parameter set k (1..2), same indexing as julia
#   DMFT_NCYCLES=n    QMC cycles per iteration (default 2000000)
#   DMFT_U=x          override U (U = 0 gives the exact cross-check)

import json
import os
import time
import numpy as np
from triqs.gf import Gf, MeshImFreq, make_gf_from_fourier
from triqs.operators import n
from triqs_cthyb import Solver

here = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(here, "data")
set_idx = int(os.environ.get("DMFT_SET", "0"))
n_cycles = int(os.environ.get("DMFT_NCYCLES", "2000000"))
seed = 777

sets = [(1.0, 1.0), (1.0, 5.0)]
run_sets = sets if set_idx == 0 else [sets[set_idx - 1]]

for (beta, U) in run_sets:
    # DMFT_U=... overrides U (including U = 0 for the exact cross-check);
    # the file tag reflects the *actual* U so that both solvers agree
    u_env = os.environ.get("DMFT_U")
    if u_env not in (None, ""):
        U = float(u_env)
    tag = f"b{int(beta)}_U{int(U)}"
    src = os.path.join(datadir, f"gtempo_nambu_{tag}.json")
    if not os.path.exists(src):
        print(f"skip β = {beta}, U = {U}: {src} not found", flush=True)
        continue
    with open(src) as f:
        gtempo = json.load(f)
    # local level ε_d = −U/2 (half filling); h_loc0 carries it directly
    eps_d = -U / 2
    n_iw_solver = 200                       # enough for tail fitting
    n_tau_solver = 6 * n_iw_solver + 2
    t2 = gtempo["t"] ** 2

    print(f"=== Nambu DMFT (cthyb): β = {beta}, U = {U}, "
          f"n_cycles = {n_cycles} ===", flush=True)

    iterations = []
    for itinfo in gtempo["iterations"]:
        it = itinfo["iter"]

        # Δ blocks from the GTEMPO export, looked up by ωₙ = Im(iωₙ);
        # frequencies beyond the GTEMPO grid use the analytic tail
        # Δuu, Δdd → t²/(iω) (anomalous blocks → 0)
        ws_all = np.array(itinfo["ws"])

        def make_lookup(key, tail):
            vals = np.array(itinfo[key])[:, 0] + 1j * np.array(itinfo[key])[:, 1]
            lut = {round(w, 9): v for w, v in zip(ws_all, vals)}
            def f(w):
                return lut.get(round(w.imag, 9), tail / w if w != 0 else 0)
            return f

        f_uu, f_dd = make_lookup("Delta_uu", t2), make_lookup("Delta_dd", t2)
        f_ud, f_du = make_lookup("Delta_ud", 0.0), make_lookup("Delta_du", 0.0)

        S = Solver(beta=beta, n_iw=n_iw_solver, n_tau=n_tau_solver,
                   gf_struct=[("nm", 2)], delta_interface=True)

        # 2×2 hybridisation on the Matsubara mesh (full grid, ±ωₙ); the four
        # GTEMPO blocks enter directly in the abstract-flavour basis
        #   c₀ = d↑ , c₁ = d↓†  (see header comment)
        iw_mesh = MeshImFreq(beta=beta, statistic="Fermion", n_iw=n_iw_solver)
        Delta_iw = Gf(mesh=iw_mesh, target_shape=(2, 2))
        for k in range(len(iw_mesh)):
            w = complex(iw_mesh[k])
            Delta_iw.data[k, 0, 0] = f_uu(w)
            Delta_iw.data[k, 1, 1] = f_dd(w)
            Delta_iw.data[k, 0, 1] = f_ud(w)
            Delta_iw.data[k, 1, 0] = f_du(w)
        S.Delta_tau["nm"] << make_gf_from_fourier(Delta_iw, n_tau=n_tau_solver)

        t0 = time.time()
        # h_loc0 = ε_d (n₀ − n₁) with ε_d = −U/2 is passed through solve (the
        # core solver owns it); U n↑n↓ = U n₀ (1 − n₁) in the Nambu basis
        S.solve(h_int=U * n("nm", 0) * (1 - n("nm", 1)),
                h_loc0=eps_d * (n("nm", 0) - n("nm", 1)),
                n_cycles=n_cycles, length_cycle=max(50, int(beta * 10)),
                n_warmup_cycles=2000, random_seed=seed)
        wall = time.time() - t0

        G_tau = S.G_tau["nm"]
        gtau = np.array(G_tau.data[:, :, :])            # (n_tau, 2, 2)
        tau = np.array([float(t) for t in G_tau.mesh])
        G_iw = S.G_iw["nm"]
        giw = np.array(G_iw.data[:, :, :])              # (n_iw, 2, 2)
        iw_pos = np.array([complex(w) for w in G_iw.mesh])

        # export the *abstract-flavour* G's; the translation to the physical
        # components (G↑↑ = G₀₀, G↓↓ = −G₁₁(β−τ), ⟨T d↑d↓⟩ = −G₀₁,
        # ⟨T d↓d↑⟩ = −G₁₀) is done in the comparison script
        iterations.append({
            "iter": it,
            "iw": [float(w.imag) for w in iw_pos],
            "G_iw_00": [[giw[k, 0, 0].real, giw[k, 0, 0].imag]
                        for k in range(len(iw_pos))],
            "G_iw_11": [[giw[k, 1, 1].real, giw[k, 1, 1].imag]
                        for k in range(len(iw_pos))],
            "G_iw_01": [[giw[k, 0, 1].real, giw[k, 0, 1].imag]
                        for k in range(len(iw_pos))],
            "G_iw_10": [[giw[k, 1, 0].real, giw[k, 1, 0].imag]
                        for k in range(len(iw_pos))],
            "G_tau_00": [float(v.real) for v in gtau[:, 0, 0]],
            "G_tau_11": [float(v.real) for v in gtau[:, 1, 1]],
            "G_tau_01": [float(v.real) for v in gtau[:, 0, 1]],
            "G_tau_10": [float(v.real) for v in gtau[:, 1, 0]],
            "tau": [float(t) for t in tau],
        })
        print(f"  iter {it:2d} done ({wall:.1f} s, "
              f"sign = {S.average_sign:.3f})", flush=True)

    out = {"solver": "cthyb", "beta": beta, "U": U, "nw": n_iw_solver,
           "n_cycles": n_cycles, "iterations": iterations}
    with open(os.path.join(datadir, f"ctqmc_nambu_{tag}.json"), "w") as f:
        json.dump(out, f)
    print(f"written ctqmc_nambu_{tag}.json", flush=True)

print("all done")
