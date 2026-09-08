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

Each script is intentionally short: the physics is set up in a few lines
(bath, model, lattice), the tensor-network parameters are collected in the
algorithm structs (`ExactTTIIF`, `TDVPIF`, truncations), and the observables
are evaluated with the functions documented in [API Reference](@ref "API Reference").
