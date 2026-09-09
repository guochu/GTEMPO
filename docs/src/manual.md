# Manual

GTEMPO solves non-equilibrium quantum impurity problems by representing the
bath hybridization influence functional (IF) as a matrix product state over
*Grassmann-valued* tensors. This page describes the building blocks and the
standard workflow; the scenario-specific walkthroughs live in the
[Practice guide](@ref "Practice guide").

## Workflow overview

```julia
bath  = fermionicbath(spectrum, β=β, μ=μ)      # 1. bath
model = AndersonIM(U=U, μ=ϵ_d)                 # 2. impurity Hamiltonian (two bands)
lat   = GrassmannLattice(N=N, δτ=δτ, contour=:imag, bands=2)   # 3. discretized contour
corr  = correlationfunction(bath, lat)         # 4. discretized hybridization
mpsI  = hybriddynamics(lat, corr, trunc=trunc) # 5. influence functional (IF)
mpsK  = sysdynamics(lat, model, trunc=trunc)   # 6. impurity dynamics (K)
mpsK  = boundarycondition!(mpsK, lat)          # 7. close the contour
obs   = gf(lat, (a, b), mpsK, mpsI)            # 8. observables
```

Every step has interchangeable implementations (steps 5 and 6 have several
algorithms, step 8 several evaluation modes), so the same script solves
Matsubara, Keldysh and Kadanoff–Baym problems by changing the `contour`
keyword of the lattice.

## Baths and correlation functions

Baths are defined by a spectral density (see `ImpurityModelBase`):
`semicircular(t)`, `spectrum(f, lb=, ub=)`, `DiracDelta(ω, α)` or a discrete
list of modes. `fermionicbath` wraps a fermionic spectrum, `bosonicbath` a
bosonic one; both accept `β` and `μ`.

`correlationfunction(bath, lattice)` discretizes the hybridization function
on the lattice, producing the coefficient matrices `η` that feed the
influence functional. On the mixed contour this yields the `τ/±`,
`τ/τ` and `±/±` blocks at once.

## Lattices

`GrassmannLattice` discretizes the contour into a chain of Grassmann
variables; each site carries an annihilation and a creation GV (`a`, `a†`).
Three contours are available:

| `contour` | branches | lattice types |
|---|---|---|
| `:imag` | `τ` | `ImagGrassmannLattice1Order/2Order` |
| `:real` / `:Keldysh` | `+`, `-` | `RealGrassmannLattice1Order/2Order` |
| `:mixed` / `:Kadanoff` | `τ`, `+`, `-` | `MixedGrassmannLattice1Order/2Order` |

Multi-band problems take a `bands` keyword (one band per orbital/spin
species). The chain arrangement of the GVs is controlled by a
*Grassmann ordering* (`ordering` keyword); see
[Grassmann lattices](@ref "Grassmann numbers and lattices (concepts)").

## Influence functional algorithms

`hybriddynamics(lattice, corr, alg)` builds the IF for four interchangeable
algorithms:

| algorithm | idea | typical use |
|---|---|---|
| `PartialIF` | product of bond-dimension-2 partial MPOs | default, robust |
| `XTRGIF` | translationally-invariant IF grown by XTRG steps | long chains |
| `ExactTTIIF` | exact exponentiation of the TI kernel | most accurate TI form |
| `TDVPIF` | second-order single-site TDVP flow | alternative to XTRGIF |

All of them take the exponential expansion algorithm (`algexpan`, an
`ExponentialExpansionAlgorithm` from ExpExp) and an MPO compression algorithm
(`algmult`, e.g. `SVDCompression`).

## Impurity dynamics

`sysdynamics(lattice, model, trunc=trunc)` evolves the impurity operator `K`.
Predefined models: `AndersonIM` (two bands, `H = μ(n₁+n₂) + U n₁n₂`),
`ToulouseIM` (single band, `U = 0`), `IRLM` and `KanamoriIM`. Custom models
subtype `ConstImpurityHamiltonian` (constant; also the quench model
`QuenchedImpurityHamiltonian`) or `AbstractTdImpurityHamiltonian`
(time-dependent, e.g. `TdImpurityHamiltonian`) and implement the
`fock_propagator` interface — `fock_propagator(model, branch, dt)` for
constant models, `fock_propagator(model, branch, dt, t)` on the real-time
branches for time-dependent ones — plus `fock_thermalstate(model, β)` if a
thermal initial state is needed; the number of bands is reported by
`num_bands`. On the mixed contour the τ leg of `K` builds the thermal density
matrix, so no separate thermal state is needed; on the real contour use
`systhermalstate!` (or start from `vacuumstate`).

## Observables

Green's functions: `gf` / `contour_ordered_gf` / `greater` / `lesser`
(integrate on the fly), the `cached_*` family (reuse precomputed
`environments`), and the `*_fast` family (MPO-based evaluation).
Occupations: `occupation`, `cached_occupation`. Currents: `electriccurrent`,
`electriccurrent_fast`, `cached_electriccurrent*`, `heatcurrent_fast`.
Density–density correlations: `nn`, `nn2`.

## Integration algorithms

`integrate(lattice, x...)` contracts GMPSs with the lattice
(`ExactIntegrate` by default, `BMPSIntegrate` for a truncated boundary-MPS
sweep). `partialintegrate` / `integrateband(s)` / `multintegrateband`
contract selected branches or bands ahead of time; these are the tools behind
the multi-orbital workflow (see
[Practice guide](@ref "Practice guide")).

## Hyperparameters and error sources

- `δτ` / `δt`: discretization of the contour; controls the discretization bias.
- `D` (`truncdimcutoff`): maximum bond dimension; the dominant controllable error.
- `ϵ`: discarded-weight cutoff of the SVD compression.
- `n`, `tol` (`OverDeterminedProny`): quality of the exponential expansion of
  the hybridization function.

Convergence should be checked by refining each of these in turn; see the
[Practice guide](@ref "Practice guide").

## Code structure

```
src/
├── lattices/            # Grassmann lattices, ContourIndex, orderings
├── grassmannmps/        # GMPS/GTerm, canonicalization, multiplication
├── mpo/                 # MPO Hamiltonians (Schur form), long-range terms
├── integration/         # lattice–GMPS contraction, cached environments
├── partialintegrate/    # partial/band integration, multi-GMPS products
├── influencefunctional/ # PartialIF / XTRGIF / ExactTTIIF / TDVPIF
├── sysdynamics/         # impurity models and their evolution
├── observables/         # GF, occupations, currents, nn
├── gvconnections/       # boundary/bulk connections
├── bcsinfluencefunctional/  # BCS bath IF
├── electronphonon/      # Fock-lattice (bosonic) path
└── auxiliary/           # misc. utilities
```
