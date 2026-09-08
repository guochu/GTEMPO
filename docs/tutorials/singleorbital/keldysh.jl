# Run from the project root with:
#   julia --project=. docs/tutorials/singleorbital/keldysh.jl

using JSON
using Serialization

using GTEMPO


# real-time (Keldysh contour) evolution up to time t, starting from the β-thermal state
function main(t; β=Inf, δt=0.1, chi=60, U=1, ϵ_d=-U/2)
	N = round(Int, t / δt)

	# fermionic bath with a semicircular spectral density
	# (β = Inf corresponds to the zero-temperature bath)
	bath = fermionicbath(semicircular(t=1), β=β, μ=0)
	# the impurity model: Anderson impurity with on-site energy ϵ_d and interaction U
	exact_model = AndersonIM(μ = ϵ_d, U=U)

	# two bands: one band per spin direction
	bands = 2
	# discretize the real-time (Keldysh) contour into N steps of size δt
	lattice = GrassmannLattice(N=N, δt=δt, contour=:Keldysh, order=1, bands=bands)
	println("number of sites, ", length(lattice))


	# discretize the bath hybridization function on the lattice
	corr = correlationfunction(bath, lattice)

	# truncation scheme: max bond dimension chi, discarded weight cutoff ϵ
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice1 = similar(lattice, bands=1)
	# @time mpsI = hybriddynamics(lattice1, corr, trunc=trunc, band=1)

	# exponential expansion of the hybridization (Prony-type fit) ...
	algexpan = OverDeterminedProny(n=20, tol=1.0e-5, verbosity=4)
	# ... feeding the translationally-invariant IF algorithm ExactTTIIF,
	# whose MPO products are compressed with SVD
	alg = ExactTTIIF(algmult=SVDCompression(trunc),algexpan=algexpan)

	# build the influence functional (IF) on the single-band lattice
	@time mpsI = hybriddynamics(lattice1, corr, alg, band=1)


	# fill (copy) the single-band IF onto each spin band
	mpsI1 = fillband(lattice, mpsI, band=1)
	mpsI2 = fillband(lattice, mpsI, band=2)

	# impurity evolution (Keldysh operator K) for the Anderson model
	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

	# close the Keldysh contour with the boundary condition per band
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end

	# precompute left/right environments shared by all observables
	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	# greater Green's function G^>(t) and lesser Green's function G^<(t)
	# for all time steps at once (fast evaluation path)
	gt = cached_greater_fast(lattice, mpsK, mpsI1, mpsI2, cache=cache)
	lt = cached_lesser_fast(lattice, mpsK, mpsI1, mpsI2, cache=cache)

	ts = [i*δt for i in 1:N]

	results = Dict("ts"=>ts, "gt"=>gt, "lt"=>lt)
	# open(data_path, "w") do f
	# 	write(f, JSON.json(results))
	# end

	return ts, gt, lt
end
