# Grassmann numbers and lattices (concepts)

This page explains the physical and combinatorial ideas behind the Grassmann
lattice; the complete function reference is in [API Reference](@ref "API Reference").
Detailed tensor conventions are described in
[Tensor and Grassmann conventions](@ref "Tensor and Grassmann conventions").

## Grassmann variables on a contour

Every discretized point of the time contour carries two Grassmann variables
(GVs), `a` and `a†`. A GV is addressed by a [`ContourIndex`](@ref):

- `j`: the time step. `j ≥ 1` are bulk time steps; `j = 0` is the traced
  boundary connecting the two ends of the contour.
- `branch`: `:τ` on the imaginary axis, `:+` / `:-` on the forward/backward
  real-time branches.
- `band`: multi-orbital problems use several independent bands.
- `conj`: `false` for `a`, `true` for `a†`.

Contour indices are totally ordered (`ContourIndex <: AbstractLatticeIndex`
implements `isless`), which is what makes contour-ordered Green's functions
well defined.

## Lattice = GV set + ordering

A `GrassmannLattice` fixes (i) the set of GVs, (ii) a *linear ordering* of
them (a `GrassmannOrdering`), and (iii) which sites are conjugation-adjacent
(`ConjugationStyle`) and how times/bands/branches interleave
(`LayoutStyle`). The ordering matters: algorithms such as
`PartialIF`/`TDVPIF` are most efficient for adjacent-conjugation orderings,
and `makestep` requires time-local orderings. Any ordering can be converted
into any other with `changeordering`, and the conversion is exact (a
permutation of site tensors, `matchindices`).

## Bands

Multi-band lattices stack `bands` copies of the single-band site pattern.
An MPS built on a one-band lattice can be inserted into a given band of a
multi-band lattice with `fillband`; bands can be permuted in place with
`swapband!` (used e.g. to build the multi-band influence functional from
band-wise two-point correlation functions).

## Influence functional and impurity dynamics

The influence functional is a Grassmann MPS `I` = exp of the bath-induced
quartic Grassmann action, discretized on the lattice. Together with the
impurity propagator `K` (a product of Fock-space propagators, one per time
step) and the boundary connection that traces the Grassmann numbers, any
correlator is a plain contraction:

```
⟨T 𝒪₁(t₁) 𝒪₂(t₂) ...⟩ = integrate(lattice, 𝒪₁𝒪₂..., I, K) / Z
```

where `Z = integrate(lattice, I, K)`. Precomputing `environments(lattice,
K, I)` turns every subsequent observable evaluation into a cheap tensor
insertion (`cached_gf`, `cached_occupation`, ...).
