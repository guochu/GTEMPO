# API Reference

## Contours and lattices

GTEMPO discretizes a time contour into a *Grassmann lattice*. Three contour
families are available:

| Contour                 | Keyword                          | Lattice type              |
|-------------------------|----------------------------------|---------------------------|
| imaginary (Matsubara)   | `contour = :imag`                | `ImagGrassmannLattice`    |
| real (Keldysh)          | `contour = :real` / `:Keldysh`   | `RealGrassmannLattice`    |
| mixed (Kadanoff–Baym)   | `contour = :mixed` / `:Kadanoff` | `MixedGrassmannLattice`   |

Each lattice is fully determined by a set of *Grassmann variables* (GVs), one
per site, addressed through [`GrassmannContourIndex`](@ref): a time step `1 ≤ j ≤ k`
on the discretized contour, a `branch` (`:τ`, `:+`, `:-`), a `band` index and
a conjugation flag distinguishing `a` from `a†`. The contour closure is not a
site: [`boundarycondition`](@ref) connects the physical end GVs directly.

```@docs
GrassmannLattice
GrassmannContourIndex
branches
indexmappings
matchindices
swapband
swapband!
fillband
vacuumstate
makestep
timesteps
```

## Orderings

The same lattice can be walked in different *Grassmann orderings*
(`GrassmannOrdering` subtypes such as `A1Ā1B1B̄1a1ā1b1b̄1`). The ordering
determines both the conjugation structure of the sites and the efficiency of
downstream algorithms:

- `ConjugationStyle`: adjacent-conjugation (`AdjacentConjugation`) versus
  general orderings (`GeneralConjugation`).
- `LayoutStyle`: how time steps, branches and bands interleave
  (`TimeLocalLayout`, `GeneralLayout`).

Orderings can be converted with `changeordering` / `toadjacentordering`,
bands can be permuted with `swapband!` and duplicated with `fillband`.

## Influence functionals

```@docs
InfluenceFunctionalAlgorithm
PartialIF
XTRGIF
ExactTTIIF
TDVPIF
hybriddynamics
influenceoperators
hybriddynamicsstepper!
```

## Impurity dynamics

```@docs
sysdynamics
sysdynamicsstepper!
sysdynamics_forward!
sysdynamics_backward!
sysdynamics_imaginary!
systhermalstate!
sysinitialstate
AndersonIM
ToulouseIM
KanamoriIM
IRLM
ImpurityHamiltonian
QuenchedImpurityHamiltonian
TdImpurityHamiltonian
TdImpurityOp
AbstractImpurityHamiltonian
ConstImpurityHamiltonian
AbstractTdImpurityHamiltonian
num_bands
```

## Boundary and bulk connections

```@docs
boundarycondition
boundarycondition!
bulkconnection
bulkconnection!
```

## Observables

### Plain evaluation

```@docs
integrate
gf
contour_ordered_gf
greater
lesser
occupation
electriccurrent
electriccurrent_fast
heatcorrelationfunction
heatcurrent_fast
nn
nn2
insert_n
```

### Cached evaluation

```@docs
environments
expectationvalue
Zvalue
cached_gf
cached_contour_ordered_gf
cached_greater
cached_lesser
cached_occupation
cached_electriccurrent
cached_electriccurrent_fast
cached_heatcurrent_fast
cached_gf_fast
cached_greater_fast
cached_lesser_fast
```

### Fast (MPO-based) observables

`electriccurrent_fast` and `cached_electriccurrent_fast` build the current
operator as a bond-dimension-2 MPO (`build_current_mpo`) instead of summing
the two-point correlators term by term.

## Entanglement entropy of a GMPS

```@docs
entanglement_entropy
entanglement_spectrum
schmidt_values
```

## Integration of GMPS products

```@docs
IntegrationAlgorithm
ExactIntegrate
BMPSIntegrate
integrateband
integratebands
multintegrateband
partialintegrate
```

## Exponential expansion algorithms

The exponential expansion of bath correlation functions is provided by the
companion package [ExpExp](https://github.com/example/ExpExp.jl) and is
re-exported here. Available algorithms: `OverDeterminedProny`,
`DeterminedProny`, `MatrixPencil`, `LeastSquareProny`; the fitting entry
point is `exponential_expansion(f, alg)`. `GenericDecayTerm` objects are
converted into lists of `ExponentialDecayTerm` by `expand_decayterm`.

## Electron–phonon path

```@docs
FockLattice
FockMPS
FockMPO
FockMatrix
fock_propagator
```
