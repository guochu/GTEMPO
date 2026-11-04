# Tutorials

Detailed scenario-based walkthroughs (normal bath, BCS bath, electron–phonon,
multi-orbital; each in imaginary, real and mixed time) are collected in the
[Practice guide](@ref "Practice guide"). Working, self-contained tutorials live in
`docs/tutorials/`. Run them with

```julia
julia --project=. docs/tutorials/singleorbital/Matsubara.jl
```

## Single orbital

| File | Contour | Content |
|---|---|---|
| `singleorbital/Matsubara.jl` | imaginary time | Matsubara Green's function of the Anderson impurity model with `ExactTTIIF`; compares against the analytical Matsubara solution |
| `singleorbital/keldysh.jl` | real time (Keldysh) | non-equilibrium real-time evolution, greater/lesser Green's functions |
| `singleorbital/kadanoff.jl` | mixed time (Kadanoff–Baym) | two-time dynamics on the Kadanoff contour |

## Multi flavor

| File | Contour | Content |
|---|---|---|
| `multiflavor/keldysh.jl` | real time | two-orbital (multi-flavor) impurity on the Keldysh contour |
| `multiflavor/kadanoff.jl` | mixed time | two-orbital impurity on the Kadanoff–Baym contour |

The `multiflavor/result/` directory contains precomputed reference data
(JSON) produced by these scripts.

## First-order IF error scaling

`firstorderif/` contains convergence benchmarks for the first-order
influence functional: for each impurity model (Toulouse `U=0`, Anderson
`U=1`) and each bath (normal `DiracDelta` single mode, BCS single mode with
gap Δ) and each contour, a self-contained script sweeps the discretization
and checks that the error against an exact-diagonalization reference obeys
first-order scaling. Each run writes `<axis>.csv` (data) and `<axis>.png`
(log–log scaling plot) next to itself; `common.jl` in each model folder
holds the shared drivers (sweep, order fit, plotting).

| Folder | Model | Bath | Scripts |
|---|---|---|---|
| `firstorderif/toulouse/` | `ToulouseIM` (U=0) | normal / BCS | `realtime.jl`, `imagtime.jl`, `mixedtime.jl`, `bcs_realtime.jl`, `bcs_imagtime.jl`, `bcs_mixedtime.jl` |
| `firstorderif/singlemodebath/` | `AndersonIM` (U=1) | normal / BCS | same six scripts |

Run one with

```julia
julia --project=. docs/tutorials/firstorderif/toulouse/realtime.jl
```

Each script prints the error table (with per-doubling orders), saves the
data, and renders the log–log plot via a bundled matplotlib snippet
(`python3` required). Measured orders are recorded in the header comment of
each script; the BCS variants also document the ↓↓-channel assembly fix
they exposed (see the repository changelog).

Each script is intentionally short: the physics is set up in a few lines
(bath, model, lattice), the tensor-network parameters are collected in the
algorithm structs (`ExactTTIIF`, `TDVPIF`, truncations), and the observables
are evaluated with the functions documented in [API Reference](@ref "API Reference").
