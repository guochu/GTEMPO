# Practice guide

Scenario-based walkthroughs with runnable code. All snippets assume
`using GTEMPO` and use small grids so they run quickly; refine `δτ`/`δt`,
`N` and the truncation `D` for production.

## 1. Normal fermionic bath

A bath defined by a spectral density (`semicircular`, `spectrum`, ...) coupled
to an Anderson impurity. This is the standard TEMPO setting.

### 1.1 Imaginary time (Matsubara)

Equilibrium problem: the τ leg of the contour builds the thermal state, and
the Matsubara Green's function `G(τ)` and partition function `Z` are the
observables.

```julia
beta = 1.0; dtau = 0.1; N = round(Int, beta/dtau)
trunc = truncdimcutoff(D=120, ϵ=1e-6)

bath  = fermionicbath(semicircular(t=1), β=beta, μ=0)
model = AndersonIM(U=1.0, μ=-0.5)

lat  = GrassmannLattice(N=N, δτ=dtau, contour=:imag)
corr = correlationfunction(bath, lat)

mpsI = hybriddynamics(lat, corr, trunc=trunc)
mpsK = sysdynamics(lat, model, trunc=trunc)
mpsK = boundarycondition!(mpsK, lat)                 # τ-periodic boundary

cache = environments(lat, mpsK, mpsI)
println("Z = ", Zvalue(cache))
gtau = cached_gf_fast(lat, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
```

The influence functional can instead be merged directly into the dynamics
state. `hybriddynamics!(mpsK, lat, corr, alg)` multiplies the IF into `mpsK`
in place (supported by `PartialIF`, `XTRGIF`, `ExactTTIIF`, `TDVPIF`), so the
observables take a single GMPS:

```julia
mpsK = sysdynamics(lat, model, trunc=trunc)
mpsK = boundarycondition!(mpsK, lat)
hybriddynamics!(mpsK, lat, corr, ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0))
Z  = integrate(lat, mpsK)
gtau = [gf(lat, (ContourIndex(i, conj=false, branch=:τ, band=1),
                 ContourIndex(1, conj=true,  branch=:τ, band=1)), mpsK, Z=Z)
        for i in 1:lat.k]
```

### 1.2 Real time (Keldysh)

Real-time evolution from a thermal initial state: greater and lesser GFs.

```julia
lat  = GrassmannLattice(N=8, δt=0.05, contour=:real)
corr = correlationfunction(bath, lat)

mpsI = hybriddynamics(lat, corr, trunc=trunc)          # or ExactTTIIF/XTRGIF/TDVPIF
mpsK = sysdynamics(lat, model, trunc=trunc)
mpsK = boundarycondition!(mpsK, lat)
mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=beta)

cache = environments(lat, mpsK, mpsI)
gt = [-im * cached_greater(lat, k, mpsK, mpsI; cache=cache) for k in 1:lat.k]
lt = [ im * cached_lesser(lat, k, mpsK, mpsI; cache=cache) for k in 1:lat.k]
```

### 1.3 Mixed time (Kadanoff–Baym)

Imaginary equilibration followed by real-time evolution: `gt`, `lt` and the
Matsubara component from one calculation.

```julia
lat  = GrassmannLattice(Nt=6, δt=0.05, Nτ=10, δτ=0.1, contour=:mixed)
corr = correlationfunction(bath, lat)
mpsI = hybriddynamics(lat, corr, trunc=trunc)
mpsK = sysdynamics(lat, model, trunc=trunc)
mpsK = boundarycondition!(mpsK, lat)

cache = environments(lat, mpsK, mpsI)
gt = [cached_greater(lat, k, mpsK, mpsI; cache=cache) for k in 1:lat.kt]
lt = [cached_lesser(lat, k, mpsK, mpsI; cache=cache)  for k in 1:lat.kt]
gτ = cached_gf_fast(lat, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
gτ[end] = 1 - gτ[1]
```

### 1.4 Currents and heat transport

```julia
j  = electriccurrent(lat, corr, mpsK, mpsI)                    # whole sequence
jc = cached_electriccurrent(lat, corr, mpsK, mpsI; cache=cache)
jf = electriccurrent_fast(lat, corr, Nt+1, mpsK, mpsI)         # MPO-based

jh = heatcurrent_fast(lat, bath, Nt+1, mpsK, mpsI)             # ω-weighted
```

Two leads (boundary-driven): sum the lead correlation functions for the IF,
evaluate each lead current with its own `corr`:

```julia
lcorr = correlationfunction(leftbath, lat)
rcorr = correlationfunction(rightbath, lat)
mpsI  = hybriddynamics(lat, lcorr + rcorr, trunc=trunc)
mpsI  = boundarycondition!(mpsI, lat)
jl = electriccurrent(lat, lcorr, mpsK, mpsI)
jr = electriccurrent(lat, rcorr, mpsK, mpsI)
```

### 1.5 Stepping (online) evolution

Instead of building the full lattice, grow it step by step:

```julia
lattice = GrassmannLattice(N=Nt, δt=dt, contour=:real, order=1)
lattice = similar(lattice, N=0)
mpsI = vacuumstate(lattice); mpsK = vacuumstate(lattice)
for k in 2:Nt+1
    lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
    mpsI  = hybriddynamicsstepper!(mpsI, lattice, corr, trunc=trunc)
    mpsI′ = boundarycondition(mpsI, lattice)
    mpsK  = sysdynamicsstepper!(mpsK, lattice, model, trunc=trunc)
    n_k   = real(occupation(lattice, k-1, mpsK, mpsI′, branch=:+))
end
```

On a second-order lattice the finalized IF must be built on a *copy*:

```julia
mpsI2 = hybriddynamicsstepper!(copy(mpsI), lattice, corr, finalize=true, trunc=trunc)
mpsI2 = boundarycondition(mpsI2, lattice)
mpsI  = hybriddynamicsstepper!(mpsI, lattice, corr, finalize=false, trunc=trunc)
```

### 1.6 Quench and time-dependent impurities

`QuenchedImpurityHamiltonian(hτ, ht)` evolves the impurity with `hτ` on the
imaginary-time branch (and builds the thermal state from it) and with `ht`
on the real-time branches:

```julia
model = QuenchedImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],   # τ leg:  μ0
                                    [tunneling(1, 1, coeff=μ1)])   # real:   μ1
```

`TdImpurityHamiltonian(hτ, ht, htt)` supports explicit time dependence on
the real-time branches: `htt` is a list of `TdImpurityOp` terms whose
coefficients are arbitrary functions of time, so the real branches evolve
with `ht + Σ op(t)` while the τ leg keeps the constant `hτ`:

```julia
model = TdImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],           # τ leg:  μ0
                              [tunneling(1, 1, coeff=μ1)],           # real:   μ1
                              [TdImpurityOp([tunneling(1, 1)],       # real:  + A·sin(ωt)
                                            t -> A * sin(ω * t))])
```

Both work with `sysdynamics` on all three contours (`sysdynamics_fast`
requires constant models). On the real contour the initial state is the
thermal state of `hτ`; on the mixed contour the τ leg builds it
automatically. Terms are set through the constructor only — there is no
`push!`.

## 2. BCS bath

A superconducting (BCS) bath couples through an anomalous pair term; the
lattice carries two bands (the impurity's Nambu pair) and the IF is built
with `orbital=1`. `bcsbath(normalbath, Δ)` wraps a normal bath with a gap `Δ`
(complex `Δ` supported).

### 2.1 Imaginary time

```julia
lat  = GrassmannLattice(N=10, δτ=0.1, contour=:imag, bands=2)
model = AndersonIM(U=1.0, μ=-ϵ_d)
mpsK = sysdynamics(lat, model, trunc=trunc)
for band in 1:2
    mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
end

normal = fermionicbath(DiracDelta(ω=1.0, α=0.5), β=β)
bath2  = bcsbath(normal, Δ=0.6)
corr   = correlationfunction(bath2, lat)
mpsI   = hybriddynamics(lat, corr, orbital=1, trunc=trunc)

cache = environments(lat, mpsK, mpsI)
gτ = cached_gf_fast(lat, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
```

Sanity check: `bcsbath(bath, Δ=0)` must reproduce two independent normal
bands — build the IF once per band (`hybriddynamics(..., band=1)` +
`swapband`) and compare.

### 2.2 Real time

Identical structure with `contour=:real`; add
`systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)` for a thermal start
and use `cached_greater_fast` / `cached_lesser_fast` for the observables.

### 2.3 Mixed time

Same workflow on `contour=:Kadanoff` (`Nt`/`δt` plus `Nτ`/`δτ`); the mixed
greater/lesser functions follow from `cached_greater` / `cached_lesser`
exactly as in the normal-bath case.

Useful checks in every contour:

* `hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc2)` (gate-by-gate
  construction) must agree with the fast construction:
  `distance(mpsI, mpsI′)/norm(mpsI) < 1e-5`.
* `Δ = 0` must reduce the problem to two independent normal bands.

## 3. Electron–phonon coupling

A bosonic (phonon) bath hybridizes on a **Fock lattice**; the fermionic
problem is carried by an ordinary Grassmann lattice, and the two are glued by
`reweighting!`.

### 3.1 Imaginary time

```julia
lattice = GrassmannLattice(N=10, δτ=0.1, contour=:imag, bands=bands)
flat    = FockLattice(N=10, δτ=0.1, contour=:imag, order=1, bands=bands)

pbath = bosonicbath(DiracDelta(ω=1.0, α=0.5), β=β)
corr  = correlationfunction(pbath, flat)
mpsI  = hybriddynamics(flat, corr, trunc=trunc)

model = (U == 0 ? ToulouseIM(μ=μ) : AndersonIM(U=U, μ=μ))  # AndersonIM is always two-band
mpsK  = sysdynamics(lattice, model, trunc=trunc)
adt   = reweighting!(lattice, mpsK, flat, mpsI, trunc=trunc)
for band in 1:bands
    adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
end

cache = environments(lattice, adt)
gτ = cached_gf_fast(lattice, adt; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
```

### 3.2 Real time

`FockLattice(..., contour=:real)` and the same `reweighting!` glue; the
observables are `cached_greater` / `cached_lesser` on the reweighted state.

### 3.3 Mixed time

`FockLattice(..., contour=:mixed)`; identical workflow.

A phonon-only alternative builds the retarded interaction directly on the
Grassmann lattice with `retardedinteractdynamics_naive` (no Fock lattice, no
reweighting); it is limited to diagonal density-type couplings.

## 4. Multi-orbital impurities

For `norb` orbitals the impurity carries `2*norb` bands (orbital × spin), and
each bath mode couples to one band. The key observation: **the IF of a
diagonal bath is the same on every band**, so it is built once on a
single-band lattice and replicated:

```julia
bath  = fermionicbath(DiracDelta(ω=1.0, α=0.5), β=β)
model = KanamoriIM(U=1.0, J=0.2, norb=2, μ=μ)          # 4 bands

lat = GrassmannLattice(Nt=Nt, δt=dt, contour=:real, bands=4)
lattice1 = similar(lat, bands=1)
corr = correlationfunction(bath, lattice1)
mpsI = hybriddynamics(lattice1, corr, trunc=trunc)
Is = [fillband(lat, mpsI, band=b) for b in 1:4]        # replicate

mpsK = sysdynamics(lat, model, trunc=trunc)
for band in 1:4
    mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
end
mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)
```

### 4.1 Observables

Direct evaluation builds a 4-band environment:

```julia
cache = environments(lat, mpsK, Is...)
gt = [-im * cached_greater(lat, k, mpsK, Is..., band=1, cache=cache) for k in 1:lat.k]
```

### 4.2 Efficient variant: band-by-band contraction

The 4-band environment is expensive. Instead, absorb the boundary condition
and one band IF into `K` at a time (`multintegrateband` multiplies two GMPSs
while integrating out a band), shrinking the lattice until a single band
remains:

```julia
algmult = DMRG1(trunc=trunc2)                # DMRG-style compression
mps_adt = mpsK
lattice_tmp = lat
for band in 1:lat.bands-1
    mps_adt = boundarycondition!(mps_adt, lattice_tmp, band=1)
    mpsI1 = fillband(lattice_tmp, mpsI, band=1)
    mps_adt = multintegrateband(lattice_tmp, mps_adt, mpsI1, algmult, band=1)
    lattice_tmp = similar(lattice_tmp, bands=lattice_tmp.bands-1)
end
mps_adt = boundarycondition!(mps_adt, lattice1, band=1)

cache = environments(lattice1, mps_adt, mpsI)
gt = cached_greater_fast(lattice1, mps_adt, mpsI, cache=cache)
lt = cached_lesser_fast(lattice1, mps_adt, mpsI, cache=cache)
```

This band-reduction trick (see `docs/tutorials/multiflavor/*.jl`, `main2`)
is the recommended way to handle multi-orbital problems.

## 5. Choosing the IF algorithm

| algorithm | construction cost | best for |
|---|---|---|
| `PartialIF` | cheap, streaming | default; short-to-medium chains |
| `XTRGIF(k=5)` | moderate | long chains, translational setting |
| `ExactTTIIF` | moderate | highest accuracy in the TI sector |
| `TDVPIF(δ=0.1)` | moderate | alternative second-order flow |

All accept `algexpan` (ExpExp algorithm) and `algmult` (MPO compression);
`verbosity > 0` prints expansion diagnostics.

## 6. Convergence rules of thumb

1. `D`: increase until observables stop moving; monitor with
   `bond_dimension(mpsI)`.
2. `δt`/`δτ`: halve and compare; second-order lattices are more accurate per
   step but cost roughly twice per step.
3. Expansion quality: `OverDeterminedProny(n=..., tol=...)` — increase `n`
   until the IF stops changing.
4. Cross-check algorithms: `PartialIF` vs `ExactTTIIF` should agree within
   truncation error (see `test/normalbath`).
