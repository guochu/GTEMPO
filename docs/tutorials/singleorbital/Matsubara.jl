# Run from the project root with:
#   julia --project=. docs/tutorials/singleorbital/Matsubara.jl

using JSON
using Serialization

using GTEMPO


function main(;β=10, δτ=0.1, chi=60,  U=1, ϵ_d=-U/2)
	N = round(Int, β / δτ)

	# fermionic bath with a semicircular spectral density, at inverse temperature β
	# and chemical potential μ (of the bath modes)
	bath = fermionicbath(semicircular(t=1), β=β, μ=0)
	# the impurity model: Anderson impurity with on-site energy ϵ_d and interaction U
	exact_model = AndersonIM(μ = ϵ_d, U=U)

	# two bands: the Anderson impurity carries one band per spin direction
	bands = 2
	# discretize the imaginary (Matsubara) contour into N steps of size δτ
	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=bands)
	println("number of sites, ", length(lattice))


	# discretize the bath hybridization function on the lattice
	corr = correlationfunction(bath, lattice)

	# truncation scheme: max bond dimension chi, discarded weight cutoff ϵ
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	# the bath hybridization is band-independent: build the influence
	# functional (IF) on a single-band lattice ...
	lattice1 = similar(lattice, bands=1)
	@time mpsI = hybriddynamics(lattice1, corr, trunc=trunc, band=1)

	# algexpan = OverDeterminedProny(n=20, tol=1.0e-5, verbosity=4)
	# alg = ExactTTIIF(algmult=SVDCompression(trunc),algexpan=algexpan)
	# @time mpsI = hybriddynamics(lattice1, corr, algmult, band=1)


	# ... and fill (copy) it onto each band of the two-band lattice
	mpsI1 = fillband(lattice, mpsI, band=1)
	mpsI2 = fillband(lattice, mpsI, band=2)

	# impurity evolution (Keldysh operator K) for the Anderson model
	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

	# close the contour: imaginary-time periodic boundary condition per band
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end

	# precompute left/right environments shared by all observables
	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	println("partition function is ", Zvalue(cache))

	# Matsubara Green's function G(τ) for all τ at once (fast evaluation path)
	gt = cached_gf_fast(lattice, mpsK, mpsI1, mpsI2; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)

	ts = [i*δτ for i in 1:N]

	results = Dict("ts"=>ts, "gtau"=>gt)
	# open(data_path, "w") do f
	# 	write(f, JSON.json(results))
	# end

	return ts, gt
end
