# GTEMPO.jl

**GTEMPO** (Grassmann TEMPO) solves quantum impurity problems by representing
the influence functional (IF) of the bath as a matrix product state over
*Grassmann-valued* tensors. Compared with conventional TEMPO-style algorithms,
the Grassmann formulation removes the need for a huge local Hilbert space and
enables exact, compact representations of the IF on arbitrary time-contour
orderings.

## How it works

A calculation is organized around a few building blocks:

1. **Bath** — define a fermionic or bosonic bath with a spectral density
   (`fermionicbath`, `bosonicbath`).
2. **Lattice** — discretize the time contour into a `GrassmannLattice`
   (`contour = :imag`, `:real` or `:mixed`), possibly with several bands.
3. **Correlation function** — `correlationfunction(bath, lattice)` computes
   the discrete hybridization matrices on the chosen contour.
4. **Influence functional** — `hybriddynamics` builds the Grassmann MPS `I`
   of the IF, with four interchangeable construction algorithms
   (`PartialIF`, `XTRGIF`, `ExactTTIIF`, `TDVPIF`).
5. **Impurity dynamics** — `sysdynamics` builds the impurity propagator `K`;
   `boundarycondition!` closes the contour, `systhermalstate!` prepares a
   thermal initial state.
6. **Observables** — evaluate Green's functions (`gf`, `greater`, `lesser`),
   occupations (`occupation`, `cached_occupation`) and currents
   (`electriccurrent`, `electriccurrent_fast`, ...) either directly or
   through precomputed environment caches (`environments`).

A minimal real-time example:

```julia
using GTEMPO

bath  = fermionicbath(semicircular(5.0), β = 1.0, μ = 0.0)
model = AndersonIM(U = 0.0, μ = 1.25π)
lattice = GrassmannLattice(N = 10, δt = 0.05, contour = :real)
trunc = truncdimcutoff(D = 100, ϵ = 1e-10)

corr = correlationfunction(bath, lattice)
I    = hybriddynamics(lattice, corr, trunc = trunc)
K    = sysdynamics(lattice, model, trunc = trunc)
K    = boundarycondition!(K, lattice)

cache = environments(lattice, K, I)
gt = [-im * cached_greater(lattice, k, K, I; cache = cache) for k in 1:lattice.k]
```

## Contents

```@contents
Pages = ["grassmann_lattice.md", "tensor_and_grassmann_conventions.md", "tutorials.md", "examples.md", "api.md"]
Depth = 2
```
