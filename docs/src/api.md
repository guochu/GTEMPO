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
per site, addressed through [`ContourIndex`](@ref): a time step `j` (0-based,
with the traced boundary at `j = 0`), a `branch` (`:τ`, `:+`, `:-`), a
`band` index and a conjugation flag distinguishing `a` from `a†`.

```@docs
GrassmannLattice
ContourIndex
index
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
(`GrassmannOrdering` subtypes such as `A1Ā1B1B̄1a1ā1b1b̄1`). The ordering
determines both the conjugation structure of the sites and the efficiency of
downstream algorithms:

- `ConjugationStyle`: adjacent-conjugation (`AdjacentConjugation`) versus
  general orderings (`GeneralConjugation`).
- `LayoutStyle`: how time steps, branches and bands interleave
  (`TimeLocalLayout`, `BandLocalLayout`, `BranchLocalLayout`,
  `GeneralLayout`).

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
KanamoriIM
IRLM
```

## Boundary and bulk connections

```@docs
boundarycondition
boundarycondition!
boundarycondition_branching
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
FockMatrix
fock_propagator
```
