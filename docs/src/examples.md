# Examples by scenario

Short, self-contained snippets covering the main use cases. Each snippet is
minimal; in production runs you would increase `N`, `δt`/`δτ` resolution and
the truncation parameters. Longer runnable scripts live in
[Tutorials](@ref "Tutorials").

```julia
using GTEMPO
```

## Imaginary-time (Matsubara) calculation

Equilibrium impurity problem on the imaginary contour: Matsubara Green's
function `G(τ)` and the partition function.

```julia
β, N = 1.0, 10                          # inverse temperature, # of τ steps
δτ = β / N
trunc = truncdimcutoff(D=120, ϵ=1e-6)   # bond dimension + cutoff

bath  = fermionicbath(semicircular(5.0), β=β, μ=0.0)  # continuous bath
model = AndersonIM(U=0.0, μ=1.25π)                    # impurity model

lat  = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
corr = correlationfunction(bath, lat)

mpsI = hybriddynamics(lat, corr, trunc=trunc)  # influence functional
mpsK = sysdynamics(lat, model, trunc=trunc)    # impurity evolution
mpsK = boundarycondition!(mpsK, lat)           # τ-periodic boundary

Z = integrate(lat, mpsK, mpsI)                 # partition function
gτ = [gf(lat, (ContourIndex(i, conj=false, branch=:τ, band=1),
                ContourIndex(1, conj=true,  branch=:τ, band=1)), mpsK, mpsI; Z=Z)
      for i in 1:lat.k]
```

## Real-time (Keldysh) dynamics

Non-equilibrium real-time evolution starting from a thermal state: greater
and lesser Green's functions.

```julia
β, δt, Nt = 1.0, 0.03, 6
trunc = truncdimcutoff(D=100, ϵ=1e-10)

bath  = fermionicbath(semicircular(5.0), β=β, μ=0.0)
model = AndersonIM(U=0.0, μ=1.25π)

lat  = GrassmannLattice(N=Nt, δt=δt, contour=:real)
corr = correlationfunction(bath, lat)

# any InfluenceFunctionalAlgorithm works; the default is PartialIF
mpsI = hybriddynamics(lat, corr, trunc=trunc)
mpsK = sysdynamics(lat, model, trunc=trunc)
mpsK = boundarycondition!(mpsK, lat)
mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)  # thermal initial state

cache = environments(lat, mpsK, mpsI)                        # shared environments
gt = [-im * cached_greater(lat, k, mpsK, mpsI; cache=cache) for k in 1:lat.k]
lt = [ im * cached_lesser(lat, k, mpsK, mpsI; cache=cache) for k in 1:lat.k]
```

Available IF algorithms (drop-in replacements in `hybriddynamics`):
`PartialIF` (default, simplest), `XTRGIF`, `ExactTTIIF`, `TDVPIF`, e.g.

```julia
mpsI = hybriddynamics(lat, corr, ExactTTIIF(algmult=SVDCompression(trunc)))
mpsI = hybriddynamics(lat, corr, XTRGIF(k=5, algmult=SVDCompression(trunc)))
mpsI = hybriddynamics(lat, corr, TDVPIF(trunc=trunc, δ=0.1))
```

## Kadanoff–Baym mixed contour

Real-time evolution preceded by an imaginary-time branch: greater, lesser
and Matsubara components from a single calculation.

```julia
Nτ, Nt = 15, 8
δτ, δt = 0.02, 0.02
β = Nτ * δτ

lat  = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)
corr = correlationfunction(bath, lat)
mpsI = hybriddynamics(lat, corr, trunc=trunc)
mpsK = sysdynamics(lat, model, trunc=trunc)
mpsK = boundarycondition!(mpsK, lat)

cache = environments(lat, mpsK, mpsI)
gt = [cached_greater(lat, k, mpsK, mpsI; cache=cache) for k in 1:lat.kt]
lt = [cached_lesser(lat, k, mpsK, mpsI; cache=cache)  for k in 1:lat.kt]

# Matsubara component, evaluated for all τ at once (fast path)
gτ = cached_gf_fast(lat, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
gτ[end] = 1 - gτ[1]
```

## Electric and heat currents

Currents on the Keldysh contour, with three interchangeable evaluation modes
(plain series, cached series, and an MPO-based "fast" variant).

```julia
# ... after the real-time setup above ...
j = electriccurrent(lat, corr, mpsK, mpsI)                # all time steps

cache = environments(lat, mpsK, mpsI)
jc = cached_electriccurrent(lat, corr, mpsK, mpsI; cache=cache)

j_fast = electriccurrent_fast(lat, corr, Nt + 1, mpsK, mpsI)  # MPO variant

# heat current: same machinery with an ω-weighted correlation function
jh_fast = heatcurrent_fast(lat, bath, Nt + 1, mpsK, mpsI)
```

Two leads (boundary-driven setup): the total influence functional is the sum
of the two lead correlation functions, and each lead current is evaluated
with its own correlation function.

```julia
leftbath  = fermionicbath(semicircular(5.0), β=β, μ= V/2)
rightbath = fermionicbath(semicircular(5.0), β=β, μ=-V/2)

lat   = GrassmannLattice(N=Nt, δt=δt, contour=:real)
lcorr = correlationfunction(leftbath, lat)
rcorr = correlationfunction(rightbath, lat)

mpsI = hybriddynamics(lat, lcorr + rcorr, trunc=trunc)
mpsI = boundarycondition!(mpsI, lat)
mpsK = sysdynamics(lat, AndersonIM(μ=-0.5, U=0), trunc=trunc)

jl = electriccurrent(lat, lcorr, mpsK, mpsI)   # current through the left lead
jr = electriccurrent(lat, rcorr, mpsK, mpsI)   # current through the right lead
```

## Stepping (online) evolution

Instead of building the full lattice at once, grow it step by step — useful
for long-time simulations and adaptive schemes.

First-order lattice:

```julia
lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1)
corr = correlationfunction(bath, lattice)

lattice = similar(lattice, N=0)        # start from an empty lattice
mpsI = vacuumstate(lattice)
mpsK = vacuumstate(lattice)
for k in 2:Nt+1
    lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
    mpsI  = hybriddynamicsstepper!(mpsI, lattice, corr, trunc=trunc)
    mpsI′ = boundarycondition(mpsI, lattice)   # non-mutating copy
    mpsK  = sysdynamicsstepper!(mpsK, lattice, model, trunc=trunc)
    n_k = real(occupation(lattice, k-1, mpsK, mpsI′, branch=:+))
    j_k = electriccurrent_fast(lattice, corr, k, mpsK, mpsI′)
end
```

Second-order lattice: the influence functional must be finalized on a *copy*
before the carried-over MPS receives the remaining increment (otherwise the
2k-2 row would be applied twice).

```julia
lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=2)
# ... same makestep loop, then inside the loop:
mpsI2 = hybriddynamicsstepper!(copy(mpsI), lattice, corr, finalize=true, trunc=trunc)
mpsI2 = boundarycondition(mpsI2, lattice)
mpsI  = hybriddynamicsstepper!(mpsI, lattice, corr, finalize=false, trunc=trunc)
mpsK  = sysdynamicsstepper!(mpsK, lattice, model, trunc=trunc)

cache = environments(lattice, mpsK, mpsI2)
# on 2Order lattices the cached observables are defined at the current step only
n_k = real(cached_occupation(lattice, mpsK, mpsI2; cache=cache))
j_k = cached_electriccurrent_fast(lattice, corr, mpsK, mpsI2; cache=cache)
```

## Few-mode and multi-orbital impurities

A few discrete bath modes (here one `DiracDelta` mode per band) coupled to a
multi-orbital impurity: each band gets its own copy of the single-band
influence functional.

```julia
bath  = fermionicbath(DiracDelta(ω=1.0, α=0.5), β=β)
model = KanamoriIM(U=1.0, J=0.2, norb=2, μ=-0.7)   # 4 impurity bands

lat = GrassmannLattice(N=Nt, δt=δt, contour=:real, bands=4)
lattice1 = similar(lat, bands=1)                   # bath lives on one band ...
corr = correlationfunction(bath, lattice1)
mpsI = hybriddynamics(lattice1, corr, trunc=trunc)
Is = [fillband(lat, mpsI, band=b) for b in 1:4]    # ... then is filled to all

mpsK = sysdynamics(lat, model, trunc=trunc)
for band in 1:4
    mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
end
mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)

cache = environments(lat, mpsK, Is...)
gt = [-im * cached_greater(lat, k, mpsK, Is..., band=1, cache=cache) for k in 1:lat.k]
lt = [ im * cached_lesser(lat, k, mpsK, Is..., band=1, cache=cache) for k in 1:lat.k]
```

Available impurity models: `AndersonIM`, `IRLM`, `KanamoriIM`.

## Electron–phonon coupling

Phonon (bosonic) hybridization is built on a Fock lattice and combined with
the fermionic problem through reweighting.

```julia
lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1, bands=bands)
flat = FockLattice(N=Nt, δt=δt, contour=:real, order=1, bands=bands)

pbath = bosonicbath(DiracDelta(ω=1.0, α=0.5), β=β)
corr  = correlationfunction(pbath, flat)      # bosonic hybridization
mpsI  = hybriddynamics(flat, corr, trunc=trunc)

model = AndersonIM(U=U, μ=μ)
mpsK = sysdynamics(lattice, model, trunc=trunc)
mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)
for band in 1:bands
    mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
end

adt = reweighting!(lattice, mpsK, flat, mpsI, trunc=trunc)  # reweighted FockMPS

cache = environments(lattice, adt)
gt = [-im * cached_greater(lattice, k, adt, band=1, cache=cache) for k in 1:Nt+1]
lt = [-im * cached_lesser(lattice, k, adt, band=1, cache=cache) for k in 1:Nt+1]
```

## GMPS integration utilities

`integrate` contracts GMPSs with the lattice. `Zvalue`-like utilities and
partial contractions help reuse intermediate results.

```julia
# product of several GMPS: exact zipup, step-by-step
Z = integrate(lat, mpsK, mpsI)                     # two GMPS
Z = integrate(lat, x, y, z)                        # any number of GMPS

# the product can also be formed explicitly and then integrated
Z2 = integrate(lat, x * y * z)

# truncated boundary-MPS integration
Z3 = integrate(lat, x, y, z; alg=BMPSIntegrate(trunc))

# partial integration: contract several GMPS on selected bands/branches
xyp = partialintegrate(lat, SVDCompression(trunc), x, y; branchs=(:τ,), bands=(1,))

# integrate out individual bands (shrinks the lattice)
xb = integrateband(lat, x; band=1)
x13 = integratebands(lat, x, (1, 3))

# multiply two GMPS while integrating out one band
xy = multintegrateband(lat, x, y; trunc=trunc)
```
