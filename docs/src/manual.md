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

Every algorithm also offers the in-place `hybriddynamics!(gmps, lattice, corr,
alg; band)`, which multiplies the influence functional directly into an
existing `GrassmannMPS` — e.g. the impurity dynamics obtained from
`sysdynamics` — so that the dynamics and the bath influence live in a single
state. `ExactTTIIF` in particular builds its IF term by term and multiplies
each term into `gmps` incrementally; on a multi-band lattice the terms are
built on a single-band lattice and expanded via `fillband`.

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

## Background theory: from path integrals to GTEMPO

This section explains *why* GTEMPO works, in plain terms, and maps each idea to
the paper that introduced it and to the code that implements it. It is written
for readers who know quantum mechanics but have not seen influence-functional
(IF) or tensor-network methods before. The papers are cited by their
bibliographic references (journal or arXiv IDs), so they can be looked up
independently of this package.

### Paper map

| paper | contribution | lives in the code as |
|---|---|---|
| R. Chen, "Path integral formalism for quantum open systems", Ann. Phys. **480**, 170083 (2025) | coherent-state path integrals, Feynman–Vernon IF, three contours, generating-functional tricks | `correlationfunction`, `fermionicbath`/`bosonicbath`, current observables |
| R. Chen, X. Xu & C. Guo, "Grassmann time-evolving matrix product operators for quantum impurity models", Phys. Rev. B **109**, 045140 (2024) | the GTEMPO algorithm: GMPS, partial IF, zipup; non-equilibrium | `hybriddynamics`, `sysdynamics`, `gf`/`greater`/`lesser` |
| R. Chen, X. Xu & C. Guo, "Grassmann time-evolving matrix product operators for equilibrium quantum impurity problems", New J. Phys. **26**, 013019 (2024) | imaginary-time GTEMPO, boundary conditions, GV ordering, cost scaling | `contour = :imag`, `boundarycondition!`, `gτ_series` |
| R. Chen, X. Xu & C. Guo, "Real-time impurity solver using Grassmann time-evolving matrix product operators", Phys. Rev. B **109**, 165113 (2024) | spectral function without analytic continuation (quench + equilibration) | real-time `systhermalstate!` workflow, `greater`/`lesser` |
| R. Chen & C. Guo, "Solving equilibrium quantum impurity problems on the L-shaped Kadanoff–Baym contour", Phys. Rev. B **110**, 165114 (2024) | GTEMPO on the mixed contour, only two controlled errors | `contour = :mixed` |
| C. Guo & R. Chen, "Infinite Grassmann time-evolving matrix product operator method in the steady state", Phys. Rev. B **110**, 045106 (2024) | iGTEMPO: infinite GMPS, cost independent of evolution time | ideas behind `XTRGIF` |
| C. Guo & R. Chen, "Infinite Grassmann time-evolving matrix product operator method for zero-temperature equilibrium quantum impurity problems", Phys. Rev. B **110**, 165119 (2024) | zero-temperature imaginary-time iGTEMPO, small bond dimension | `XTRGIF`, small-`χ` imaginary runs |
| Z. Sun, R. Chen, Z. Li & C. Guo, "Infinite Grassmann time-evolving matrix product operators for quantum impurity problems after a quench", Phys. Rev. B **112**, 125145 (2025) | non-equilibrium iGTEMPO: window GMPS for the quench | quenched/time-dependent models |
| Z. Sun, R. Chen, Z. Li & C. Guo, "Scalable tensor network algorithm for quantum impurity problems", Phys. Rev. B **112**, 155115 (2025) | multiflavor GTEMPO: integrate out unobserved flavors first | multi-orbital workflow, `partialintegrate` |
| Z. Sun, Z. Li & C. Guo, "Efficient and accurate tensor network algorithm for Anderson impurity problems", arXiv:2510.11459 (2025) | exact TTI IF via exponential fits of the hybridization | `ExactTTIIF` |
| C. Guo, W. Wu, X. Xu, P.-X. Chen, C. Yue, T. Jiang & R. Chen, "Grassmann time-evolving matrix product operators for superconducting quantum impurity model", arXiv:2604.23301 (2026) | Nambu–GTEMPO for BCS baths via Bogoliubov transformation | `bcsbath`, `bcsinfluencefunctional/` |
| R. Chen, L. Gu & C. Guo, "Tensor network algorithm to solve polaron impurity problems", Chin. Phys. Lett. **42**, 120701 (2025) | electron + phonon baths both integrated out | `bosonicbath`, `FockLattice`, `reweighting!` |

### 1. Where the IF comes from (Chen 2025)

A quantum open system is "impurity + bath". The bath is large and
non-interacting, so it can be *integrated out analytically*: the Feynman–Vernon
influence functional (IF) `I` is the only memory the bath leaves on the
impurity. GTEMPO never stores the bath wave function — the bath enters the
calculation only through a *hybridization function* `Δ(τ,τ')` fixed by its
spectral density.

For **bosons** the path integral is over ordinary complex numbers (coherent
states are the eigenstates of `b̂`). For **fermions** the coherent-state
eigenvalues must anticommute, so they are **Grassmann numbers**: `θθ' = -θ'θ`,
`θ² = 0`. This single difference is why fermionic impurity problems resisted
TEMPO for so long — and why the code in this package is built on
`Z2Tensors.jl`, a tensor library whose indices carry a fermionic (Z₂) parity.

On a discretized contour the IF factorizes (for a fermionic bath):

```
I[ā,a] = exp(-Σ_{j,k} ā_j η_{j,k} a_k) ,
```

where the matrix `η` is the discretized hybridization. This `η` is exactly
what `correlationfunction(bath, lattice)` produces from a spectral density:

```julia
bath = fermionicbath(semicircular(t=1), β=1.0, μ=0)   # spectral density J(ε)
lat  = GrassmannLattice(N=10, δτ=0.1, contour=:imag)
corr = correlationfunction(bath, lat)                  # η matrices on the contour
```

Three contours are supported — imaginary time (Matsubara), real time
(Keldysh) and the L-shaped Kadanoff–Baym contour (`contour = :imag`, `:real`,
`:mixed`) — corresponding to the three cases worked out in the review.

### 2. Turning the path integral into a tensor network (Chen, Xu & Guo 2024a, 2024b)

The full path integral is a product of three factors:

```
Z = ∫ D[ā,a] K[ā,a] I[ā,a]         (K: impurity dynamics, I: bath influence)
```

GTEMPO represents both factors as **Grassmann matrix product states (GMPS)** —
an MPS whose site tensors are Z₂-graded, so that the fermionic signs are
handled by the tensor algebra instead of by explicit swap gates.

- **`hybriddynamics(lat, corr, alg)`** builds the IF `I` as a GMPS. The
  default `PartialIF` implements the *partial-IF* construction of the papers:
  `I = exp(-Σ ā η a)` is regrouped as a product of O(N) *partial IFs*, each
  of bond dimension 2, which are then multiplied together.
- **`sysdynamics(lat, model, trunc)`** builds the impurity propagator `K`
  from the impurity Hamiltonian, and **`boundarycondition!`** closes the
  contour (the fermionic sign of the trace `⟨-a,0|...` is built in).
- The final expectation value is computed on the fly (the *zipup* algorithm):
  `K` and the per-flavor `I`'s are kept separate and contracted together only
  when an observable is evaluated. This avoids ever forming the large
  augmented density tensor `A = K·ΠI`.

```julia
bath  = fermionicbath(semicircular(5.0), β=1.0, μ=0.0)
model = AndersonIM(U=2.0, μ=-1.0)
lattice = GrassmannLattice(N=10, δt=0.05, contour=:real, bands=2)
trunc = truncdimcutoff(D=100, ϵ=1e-9)

corr = correlationfunction(bath, lattice)
I = hybriddynamics(lattice, corr, ExactTTIIF(algmult=SVDCompression(trunc)))
K = sysdynamics(lattice, model, trunc=trunc)
K = boundarycondition!(K, lattice)
gt = [greater(lattice, t, K, I) for t in 1:lattice.k]
```

Key practical lessons from the papers, all built into the package:

- **GV ordering matters.** The arrangement of the Grassmann variables along
  the chain controls the bond dimension. The time-local orderings used by
  default are the ones the papers found to be almost always optimal.
- **Imaginary time needs more bond dimension than real time.** The Matsubara
  hybridization contains growing as well as decaying exponentials, so `χ`
  grows roughly linearly with `β`; on the real-time axis `χ` saturates.
- **Only two sources of error**: time discretization (`δτ`, `δt`) and bond
  truncation (`D`, `ϵ`). This is the same error structure analyzed in the
  Kadanoff–Baym paper (Chen & Guo 2024g).

### 3. Time-translational invariance and the infinite-GTEMPO family (Guo & Chen 2024e, 2024f; Sun et al. 2025a)

Because the bath is time-independent, the hybridization depends only on the
*time difference*: `η_{j,k} = η_{j-k}`. The IF is **time-translationally
invariant (TTI)**, and its GMPS is periodic — one needs to store only a single
*unit cell*. The iGTEMPO papers exploit this to remove the total evolution
time `N` from the storage and construction cost:

- **Steady state** (Guo & Chen 2024e): at long times the initial state is
  forgotten, so the dynamics is translationally invariant and the IF is
  represented by an infinite GMPS. Observables are computed by first finding
  the dominant left/right eigenvectors of the single-cell transfer matrix and
  using them as boundary conditions.
- **Zero temperature** (Guo & Chen 2024f): on the imaginary axis at `β = ∞`
  every contribution of the hybridization decays exponentially (`|λ| < 1` in
  the Prony expansion), so the bond dimension stays small — the regime that is
  hardest for CTQMC.
- **After a quench** (Sun et al. 2025a): the problem has no TTI, but the
  impurity dynamics `K` is split into a *window* GMPS (the quench, finite)
  sandwiched between two infinite GMPS built with the pre-quench Hamiltonian.
  The bath influence is still translationally invariant, so the cost of the
  expensive parts stays independent of the evolution time.

The IF construction shared by these papers has two ingredients:

1. a **Prony (exponential) expansion** of the hybridization,
   `η_x ≈ Σ_l α_l λ_l^|x|` — this is `algexpan` (`OverDeterminedProny` from
   ExpExp) in the algorithm structs;
2. a **doubling trick** (analogous to XTRG): build `e^{δF}` for a tiny step
   `δ = 1/2^k` and square the infinite GMPS `k` times to reach `I = e^F`.

The `XTRGIF` algorithm in this package implements this translationally
invariant construction (it is named after the XTRG-style doubling):

```julia
I = hybriddynamics(lattice, corr, XTRGIF(algexpan=OverDeterminedProny(n=15, tol=1e-4),
                                         algmult=SVDCompression(trunc), k=5))
```

The `ExactTTIIF` algorithm is the follow-up construction of Sun, Li & Guo
(2025c): with the hybridization fitted by `n` exponentials, each exponential
term becomes a *bond-dimension-2* GMPS for which the first-order expansion is
exactly correct (a special property of Grassmann variables — the same result
fails for bosons), so the IF is built by multiplying `2n` (imaginary) or
`8n` (real-time) small GMPS with no expansion bias.

### 4. Extensions implemented in the package

- **Spectral functions without analytic continuation** (Chen, Xu & Guo 2024c):
  start from a separable initial state, evolve for an equilibration time `t₀`,
  and Fourier-transform the real-time Green's function directly. In the code
  this is the real-contour workflow with `systhermalstate!`:

  ```julia
  K = sysdynamics(lattice, model, trunc=trunc)
  K = boundarycondition!(K, lattice)
  K = systhermalstate!(K, lattice, model, trunc=trunc, β=β)   # thermal initial state
  gt = [greater(lattice, t, K, I) for t in 1:lattice.k]
  ```

- **Mixed (Kadanoff–Baym) contour** (Chen & Guo 2024g): `contour = :mixed`
  gives Matsubara and real-time Green's functions from a single calculation,
  with only discretization and truncation errors — no analytic continuation,
  no equilibration-time guess.

- **Multi-flavor (multi-orbital) impurities** (Sun, Chen, Li & Guo 2025b):
  when each flavor couples to its own bath (diagonal hybridization), the IF
  factorizes and the unobserved flavors can be integrated out *first*, cutting
  the exponential bond-dimension growth to `O(2^n χ)`. The package exposes the
  building blocks (`partialintegrate`, `integrateband(s)`, `multintegrateband`)
  used by the multi-orbital workflow in the
  [Practice guide](@ref "Practice guide").

- **Superconducting (BCS) baths** (Guo et al. 2026, Nambu–GTEMPO): a
  Bogoliubov transformation `α = u c - v c†` turns the BCS bath into a normal
  one (quasiparticle energy `ξ_k = √(ε_k² + |Δ|²)`), so the IF keeps its
  quadratic form with a 2×2 Nambu hybridization matrix. The code wraps this in
  `bcsbath`:

  ```julia
  bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=0.6)
  corr  = correlationfunction(bath2, lattice)
  I     = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
  ```

- **Electron–phonon (polaron) problems** (Chen, Gu & Guo 2025): both the
  fermionic bath and the bosonic phonon bath are integrated out. The phonon
  IF is a *retarded interaction* that is simplest on the occupation-number
  basis, so the package builds it on a `FockLattice` and reweights it onto the
  Grassmann lattice:

  ```julia
  pbath = bosonicbath(DiracDelta(ω=ω₀, α=α), β=β)   # phonon bath
  flat  = FockLattice(N=N, δτ=δτ, contour=:imag)
  mpsI  = hybriddynamics(flat, correlationfunction(pbath, flat), trunc=trunc)
  mpsK  = sysdynamics(lattice, model, trunc=trunc)
  adt   = reweighting!(lattice, mpsK, flat, mpsI, trunc=trunc)
  ```

  A subtle but important point from the polaron paper: the phonon IF
  *cannot* be obtained by naively replacing `n̂ → āa` in `e^{αn̂}`; the
  correct matrix element is `⟨a'|e^{αn̂}|a⟩ = e^{e^α ā'a}` (and not
  `e^{α ā'a}`). The Fock-lattice path handles this exactly.

### Suggested reading order

1. Chen (2025) — read Sections on the IF and the three contours; skip the
   generating-functional part on first pass.
2. Chen, Xu & Guo (2024a) — the core GTEMPO algorithm (sections on GMPS,
   partial IF and zipup).
3. Chen, Xu & Guo (2024b) — imaginary-time subtleties (boundary conditions,
   GV ordering).
4. Guo & Chen (2024e, 2024f) and Sun et al. (2025a) — the infinite-GTEMPO
   family; the papers that motivate `XTRGIF` and the long-chain regime.
5. The extension papers (2024c, 2024g, 2025b, 2025c, 2026, polaron) — read
   whichever matches your application.

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
