# reference-solution helpers shared by the model tests. `μ` is the chemical
# potential handed to AndersonIM / KanamoriIM (H_imp = μ Σ n_i).

using LinearAlgebra: exp

# ED Matsubara Green's function ⟨a(τ) a†⟩ of the full ED Hamiltonian
gτ_ed(H, a, adag, τs, β) = correlation_2op_1τ(H, a, adag, τs, β=β)

"""
    greater_lesser_ed(H, a, adag, H0, ts, β)

ED real-time greater/lesser functions with initial density matrix
`exp(-β H0)` (H0 = the decoupled Hamiltonian for a separable thermal state).
"""
function greater_lesser_ed(H, a, adag, H0, ts, β)
	ρ = exp(-β * H0)
	cache = eigencache(H)
	gt = -im .* correlation_2op_1t(H, a, adag, ρ, ts, cache, reverse=false)
	lt = im .* correlation_2op_1t(H, adag, a, ρ, ts, cache, reverse=true)
	return gt, lt
end

# Kadanoff variant: the τ leg builds the interacting thermal state exp(-β H)
function greater_lesser_ed_mixed(H, a, adag, ts, β)
	ρ = exp(-β * H)
	cache = eigencache(H)
	gt = -im .* correlation_2op_1t(H, a, adag, ρ, ts, cache, reverse=false)
	lt = im .* correlation_2op_1t(H, adag, a, ρ, ts, cache, reverse=true)
	return gt, lt
end

"""
    fermionic_setup(lattice, bath, model, trunc; β=nothing)

Standard Grassmann-path setup for a fermionic bath: build the single-band
influence functional, replicate it with `fillband`, build the impurity
dynamics, apply boundary conditions and (optionally) the thermal initial
state. Returns (mpsK, Is).
"""
function fermionic_setup(lattice, bath, model, trunc; β=nothing)
	corr = correlationfunction(bath, lattice)
	mpsI1 = hybriddynamics(lattice, corr, band=1, trunc=trunc)
	Is = [mpsI1]
	for b in 2:lattice.bands
		push!(Is, swapband(Is[1], lattice, 1, b, trunc=trunc))
	end
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	if !isnothing(β)
		mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)
	end
	return mpsK, Is
end

# ---------------------------------------------------- observable series
# All series work for any Grassmann ordering. When `cache` is given, the
# (faster) cached evaluation is used; otherwise the plain functions with an
# explicit `Z` are used (cached evaluation requires AdjacentConjugation).

# Matsubara Green's function on the imaginary-time lattice
function gτ_series(lattice, mpsK, Is...; cache=nothing)
	if isnothing(cache)
		# plain evaluation accepts a complex Z (e.g. complex BCS gaps)
		Z = integrate(lattice, mpsK, Is...)
		b = ContourIndex(1, conj=true, branch=:τ, band=1)
		return [gf(lattice, (ContourIndex(i, conj=false, branch=:τ, band=1), b), mpsK, Is..., Z=Z) for i in 1:lattice.k]
	end
	return cached_gf_fast(lattice, mpsK, Is...; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
end

# greater/lesser series on the Keldysh contour
function gtlt_series(lattice, mpsK, Is...; cache=nothing, band=1)
	if isnothing(cache)
		Z = integrate(lattice, mpsK, Is...)
		gt = [-im * gf(lattice, (ContourIndex(k, conj=false, branch=:+, band=band), ContourIndex(1, conj=true, branch=:+, band=band)), mpsK, Is..., Z=Z) for k in 1:lattice.k]
		lt = [-im * gf(lattice, (ContourIndex(1, conj=true, branch=:-, band=band), ContourIndex(k, conj=false, branch=:+, band=band)), mpsK, Is..., Z=Z) for k in 1:lattice.k]
	else
		gt = [-im * cached_greater(lattice, k, mpsK, Is..., band=band, cache=cache) for k in 1:lattice.k]
		lt = [-im * cached_lesser(lattice, k, mpsK, Is..., band=band, cache=cache) for k in 1:lattice.k]
	end
	return gt, lt
end

# greater/lesser/Matsubara series on the Kadanoff contour
function gtltgτ_series(lattice, mpsK, Is...; cache=nothing, band=1)
	if isnothing(cache)
		Z = integrate(lattice, mpsK, Is...)
		gt = [-im * gf(lattice, (ContourIndex(k, conj=false, branch=:+, band=band), ContourIndex(1, conj=true, branch=:+, band=band)), mpsK, Is..., Z=Z) for k in 1:lattice.kt]
		lt = [im * gf(lattice, (ContourIndex(1, conj=true, branch=:-, band=band), ContourIndex(k, conj=false, branch=:+, band=band)), mpsK, Is..., Z=Z) for k in 1:lattice.kt]
		gτ = [gf(lattice, (ContourIndex(k, conj=false, branch=:τ, band=band), ContourIndex(1, conj=true, branch=:τ, band=band)), mpsK, Is..., Z=Z) for k in 1:lattice.kτ]
	else
		gt = [-im * cached_greater(lattice, k, mpsK, Is..., cache=cache) for k in 1:lattice.kt]
		lt = [im * cached_lesser(lattice, k, mpsK, Is..., cache=cache) for k in 1:lattice.kt]
		gτ = cached_gf_fast(lattice, mpsK, Is...; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
		gτ[end] = 1 - gτ[1]
	end
	return gt, lt, gτ
end

# the four influence-functional construction algorithms (imaginary/real time)
if_algs(trunc) = [
	("PartialIF",  PartialIF(trunc=trunc)),
	("XTRGIF",     XTRGIF(k=3, algmult=SVDCompression(trunc), verbosity=0)),
	("ExactTTIIF", ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0)),
	("TDVPIF",     TDVPIF(trunc=trunc, δ=0.1, verbosity=0)),
]

# analytic independent-bosons reference with the μ convention (ϵ_d = μ)
independentbosons_Gτ_μ(spec; β, μ, Nτ, U=0, bands=1) =
	independentbosons_Gτ(spec, β=β, ϵ_d=μ, Nτ=Nτ, U=U, bands=bands)
independentbosons_greater_μ(spec, t; β, μ, U=0, bands=1) =
	independentbosons_greater(spec, t, β=β, ϵ_d=μ, U=U, bands=bands)
independentbosons_lesser_μ(spec, t; β, μ, U=0, bands=1) =
	independentbosons_lesser(spec, t, β=β, ϵ_d=μ, U=U, bands=bands)
