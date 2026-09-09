# Internals

This page documents the internal design of GTEMPO: the Grassmann tensor
formalism, how the lattice chain is laid out, and the index-shift convention
used by `hybriddynamics` on the imaginary-time branch — including the
numerical evidence that justifies it. See also
[Tensor conventions](@ref "Tensor and Grassmann conventions") for the
`@tensor`/`@grassmann` usage rules and
[Grassmann lattices](@ref "Grassmann numbers and lattices (concepts)") for the
combinatorial picture.

## 1. Module structure

```
src/
├── lattices/                 GrassmannLattice, ContourIndex, branch/ordering logic
├── grassmannmps/             GrassmannMPS, GTerm, canonicalization, mult!
├── mpo/mpohamiltonian/       MPO Hamiltonians (SchurMPOTensor), long-range terms
├── integration/              lattice contraction, TwosideExpectationCache
├── partialintegrate/         partialintegrate, integrateband(s), multintegrateband
├── influencefunctional/      PartialIF / XTRGIF / ExactTTIIF / TDVPIF
├── sysdynamics/              impurity models (AndersonIM, ToulouseIM, IRLM,
│                             KanamoriIM, ImpurityHamiltonian, Quenched/
│                             TdImpurityHamiltonian, ...)
├── observables/              gf / cached_* / *_fast, occupations, currents
├── gvconnections/            boundarycondition, bulkconnection
├── bcsinfluencefunctional/   BCS bath influence functionals
└── electronphonon/           Fock-lattice (bosonic) path, reweighting
```

## 2. Grassmann tensor design

### 2.1 Grassmann numbers and Z2-graded spaces

Every lattice site carries two Grassmann variables: an annihilation-like GV
`a` (`conj = false`) and a creation-like GV `a†` (`conj = true`). Tensors over
these variables are `Z2Tensors.TensorMap`s: each index is graded even/odd by
the fermion parity, stored in the `FusionBlockStructure` of the tensor. A
`GrassmannMPS` site tensor therefore has a `2 → 1` (annihilation) or `1 → 2`
(creation) flow, and the vacuum state has a definite parity on every bond.

Because the physical spaces are Z2-graded, contractions never need explicit
fermionic swap gates *within* a tensor: the block structure enforces parity
conservation, and only the *ordering of sites along the chain* matters for
signs.

### 2.2 Chain layout and orderings

A lattice is a chain of these GVs. Which physical time step / branch / band
lands on which chain position is fixed by a `GrassmannOrdering` (e.g.
`A1Ā1B1B̄1a1ā1b1b̄1` for the real-time default). The ordering determines:

* the `conj` alternation (`ConjugationStyle`: adjacent conjugations vs general),
* how time steps, branches and bands interleave (`LayoutStyle`:
  `TimeLocalLayout`, `BandLocalLayout`, `BranchLocalLayout`, `GeneralLayout`).

Along the imaginary branch the chain runs **backwards in time** (τ = β sits
next to the boundary, τ = δτ next to the real branches); real branches run
forward, the − branch runs backward in contour order.

### 2.3 Signs

The vacuum expectation value of a `GTerm` (a product of GVs at given chain
positions) is evaluated by contracting the chain from both ends into the
vacuum (`TwosideExpectationCache`). The fermionic signs of the permutation
that brings the GVs next to each other are handled internally by
`compensate_twists!` / `gpermute` during the construction of the MPO gates,
so that user-level code works with plain Grassmann products
(`GTerm(pos1, pos2, coeff=η)` → `exp(GTerm(...))`).

## 3. IF construction and the one-site τ shift

### 3.1 The convention

`hybriddynamics` for the mixed contour places each influence-functional gate
`exp(η(τ_i, τ_j) δτ² d† d)` on the chain as follows
(`partialif/mixedtime.jl`):

* the gate row of the τ branch `i` is placed at lattice position `i + 1`,
* the τ columns `j` are placed at lattice positions `j + 1`,
* real-branch (+/−) rows and columns are placed **without** a shift.

Equivalently: the correlation index `i` of every branch maps to the lattice
site `i + 1` for τ but to the site `i` for the real branches. The site `i=0`
of every branch is the boundary/junction GV glued by `boundarycondition!`
and receives no IF gate.

### 3.2 Why the shift exists

The τ branch of the contour starts at τ = 0, which is the traced junction of
the Kadanoff contour — there is no Grassmann variable there that the IF
could couple to. The discretized hybridization on the τ branch is a
*half-open* product over τ ∈ [δτ, β]:

```
I_τ = exp( Σ_{i,j=1..Nτ} η(τ_i, τ_j) δτ²  d†(τ_i) d(τ_j) )
```

Placing the gate for the coefficient `η(τ_i, τ_j)` at the chain sites of
τ_{i+1} / τ_{j+1} implements exactly this half-open product. The real
branches carry ordinary time-ordered evolution from t = 0, so their IF gates
sit on the sites of the same time index (no shift).

### 3.3 Does the shift hurt the mixed GF?

The shift is *asymmetric* between the τ branch and the real branches, so one
might worry that cross-branch Green's functions G(τ, t) come out wrong.
`performance/mixedgf/mixedgf.jl` checks this against exact diagonalization
(single- and two-band Anderson impurity, U = 0 and U ≠ 0) on the **full**
(τ_i, t_j) grid:

| branch pair | worst |GTEMPO − ED| (δτ = 0.1) | at δτ = 0.05 |
|---|---|---|
| creator on + branch | 1.8e-3 | 1.2e-3 |
| creator on − branch | 1.0e-2 | 9.8e-3 |

The deviations stay at truncation level (D = 80) and do **not** scale with
δτ — the (τ = 0 junction, t_max) corner where the − branch value is largest
is an end-point/truncation effect, not a shift artifact. Pure-τ and pure-real
GFs are exact to the same level. **The one-site τ shift is therefore correct
as implemented and requires no adjustment in `hybriddynamics`.**

### 3.4 Reading cross-branch values

Because the raw `gf` values are plain chain contractions (no re-ordering),
the physical identification of a cross-branch pair follows the contour
order τ < + < − with the calibrated conj convention (annihilator-first on
τ/+, creator-first on −):

* `gf((τ_i, false), (br, j, true))` = ⟨d(τ_i) d†(t_j)⟩  — the mixed GF above;
* `gf((br, j, false), (τ_i, true))` = −⟨d†(τ_i) d(t_j)⟩ + (anti-commutator
  support only at equal times), i.e. the lesser-type quantity.

`contour_ordered_gf` applies the standard T̂ re-ordering automatically; use it
when the ordering of the two arguments is not known a priori.

## 4. Cached evaluation and scaling

`environments(lattice, K, I...)` builds left/right chain environments with a
running normalization (`hleft_scaling` / `hright_scaling`); every expectation
value is `contract_center(...) / Zvalue(cache)` with the scalings folded in.
The `cached_gf_fast` family reuses the environments across a whole sweep of
time indices, and the `*_fast` observables replace the sum over two-point
correlators by a bond-dimension-2 MPO (`build_current_mpo`).

## 5. Sources of error

| source | controlled by | symptom |
|---|---|---|
| IF discretization | δτ, δt | O(δ) bias of correlators |
| SVD truncation | D, ϵ | loss of accuracy at large t / β |
| exponential expansion | n, tol (OverDeterminedProny) | wrong long-time decay of the IF |
| band-dimension of fast observables | fixed at 2 | exact for two-point observables |
