push!(LOAD_PATH, "../../src")


using JSON
using Serialization

using GTEMPO


function main(t; β=Inf, δt=0.1, chi=60, U=1, ϵ_d=-U/2)
	N = round(Int, t / δt)

	bath = fermionicbath(semicircular(t=1), β=β, μ=0)
	exact_model = AndersonIM(μ = ϵ_d, U=U)

	bands = 2
	lattice = GrassmannLattice(N=N, δt=δt, contour=:Keldysh, order=1, bands=bands)
	println("number of sites, ", length(lattice))


	corr = correlationfunction(bath, lattice)

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	lattice1 = similar(lattice, bands=1)
	# @time mpsI = hybriddynamics(lattice1, corr, trunc=trunc, band=1)

	algmult = ExactTranslationInvariantIF(algmult=SVDCompression(trunc=trunc))
	algexpan = PronyExpansion(n=20, tol=1.0e-5, verbosity=4)
	@time mpsI = hybriddynamics(lattice1, corr, algmult, band=1)


	mpsI1 = fillband(lattice, mpsI, band=1)
	mpsI2 = fillband(lattice, mpsI, band=2)

	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end

	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	gt = cached_greater_fast(lattice, mpsK, mpsI1, mpsI2, cache=cache)
	lt = cached_greater_fast(lattice, mpsK, mpsI1, mpsI2, cache=cache)

	ts = [i*δt for i in 1:N]

	results = Dict("ts"=>ts, "gt"=>gt, "lt"=>lt)
	# open(data_path, "w") do f
	# 	write(f, JSON.json(results))
	# end

	return ts, gt, lt
end