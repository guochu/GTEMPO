# Run from the project root with:
#   julia --project=. docs/tutorials/singleorbital/kadanoff.jl

using JSON
using Serialization

using GTEMPO


# real-time evolution up to time t after an imaginary-time equilibration
# of inverse temperature β (Kadanoff–Baym / mixed contour)
function main(t; β=10, δτ=0.1, δt=0.1, chi=60,  U=1, ϵ_d=-U/2)
	Nt = round(Int, t / δt)
	Nτ = round(Int, β/δτ)

	# fermionic bath with a semicircular spectral density
	bath = fermionicbath(semicircular(t=1), β=β, μ=0)
	# the impurity model: Anderson impurity with on-site energy ϵ_d and interaction U
	exact_model = AndersonIM(μ = ϵ_d, U=U)

	# two bands: one band per spin direction
	bands = 2
	# discretize the mixed (Kadanoff) contour: imaginary branch (Nτ, δτ)
	# followed by the real branch (Nt, δt)
	lattice = GrassmannLattice(Nt=Nt, Nτ=Nτ, δτ=δτ, δt=δt, contour=:Kadanoff, order=1, bands=bands)
	println("number of sites, ", length(lattice))


	# discretize the bath hybridization function on the lattice
	corr = correlationfunction(bath, lattice)

	# truncation scheme: max bond dimension chi, discarded weight cutoff ϵ
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	# the bath hybridization is band-independent: build the influence
	# functional (IF) on a single-band lattice ...
	lattice1 = similar(lattice, bands=1)
	@time mpsI = hybriddynamics(lattice1, corr, trunc=trunc, band=1)

	# ... and fill (copy) it onto each band of the two-band lattice
	mpsI1 = fillband(lattice, mpsI, band=1)
	mpsI2 = fillband(lattice, mpsI, band=2)

	# impurity evolution (Keldysh operator K) for the Anderson model
	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

	# close the mixed contour with the boundary condition per band
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end

	# precompute left/right environments shared by all observables
	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	# greater Green's function G^>(t) and lesser Green's function G^<(t)
	# for all real-time steps at once (fast evaluation path)
	gt = cached_greater_fast(lattice, mpsK, mpsI1, mpsI2, cache=cache)
	lt = cached_lesser_fast(lattice, mpsK, mpsI1, mpsI2, cache=cache)

	ts = [i*δt for i in 1:N]

	results = Dict("ts"=>ts, "gt"=>gt, "lt"=>lt)
	# open(data_path, "w") do f
	# 	write(f, JSON.json(results))
	# end

	return ts, gt, lt
end
