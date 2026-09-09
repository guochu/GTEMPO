@testset "API: observables" begin
	# free impurity (U = 0) coupled to a semicircular bath
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	model = ToulouseIM(μ=-0.5)
	bath = fermionicbath(spectrum_func(), β=1.0, μ=0)

	# imaginary time
	lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	corr = correlationfunction(bath, lat)
	I = hybriddynamics(lat, corr, trunc=trunc)
	K = sysdynamics(lat, model, trunc=trunc)
	K = boundarycondition!(K, lat)
	Z = integrate(lat, K, I)
	@test Z ≈ real(Z) atol=1.0e-8
	cache = environments(lat, K, I)
	@test Zvalue(cache) ≈ Z rtol=1.0e-8
	g1 = cached_gf_fast(lat, K, I; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
	@test length(g1) == lat.k
	g2 = [gf(lat, (ContourIndex(i, conj=false, branch=:τ, band=1), ContourIndex(1, conj=true, branch=:τ, band=1)), K, I; Z=Z) for i in 1:lat.k]
	@test relerr(g1, g2) < 1.0e-8
	# equilibrium occupation is time independent
	occs = [occupation(lat, i, K, I, Z=Z) for i in 1:lat.k]
	@test maximum(abs.(occs .- occs[1])) < 1.0e-8

	# real time (thermal initial state)
	lat = GrassmannLattice(N=4, δt=0.05, contour=:real)
	corr = correlationfunction(bath, lat)
	I = hybriddynamics(lat, corr, trunc=trunc)
	K = sysdynamics(lat, model, trunc=trunc)
	K = boundarycondition!(K, lat)
	K = systhermalstate!(K, lat, model, trunc=trunc, β=1.0)
	cache = environments(lat, K, I)
	g1 = [cached_gf(lat, (ContourIndex(k, conj=false, branch=:+, band=1), ContourIndex(1, conj=true, branch=:+, band=1)), K, I; cache=cache) for k in 1:lat.k]
	g2 = [gf(lat, (ContourIndex(k, conj=false, branch=:+, band=1), ContourIndex(1, conj=true, branch=:+, band=1)), K, I; Z=Zvalue(cache)) for k in 1:lat.k]
	@test relerr(g1, g2) < 1.0e-8
	gfast = cached_gf_fast(lat, K, I, c1=false, c2=true, b1=:+, b2=:+, cache=cache)
	@test length(gfast) == lat.k
	n1 = [cached_occupation(lat, i, K, I, cache=cache) for i in 1:lat.k-1]
	@test all(isreal, n1)

	# mixed time (Kadanoff)
	lat = GrassmannLattice(Nt=3, δt=0.05, Nτ=4, δτ=0.25, contour=:mixed)
	corr = correlationfunction(bath, lat)
	I = hybriddynamics(lat, corr, trunc=trunc)
	K = sysdynamics(lat, model, trunc=trunc)
	for band in 1:lat.bands
		K = boundarycondition!(K, lat, band=band)
	end
	cache = environments(lat, K, I)
	g1 = [cached_greater(lat, k, K, I, cache=cache) for k in 1:lat.kt]
	g2 = [gf(lat, (ContourIndex(k, conj=false, branch=:+, band=1), ContourIndex(1, conj=true, branch=:+, band=1)), K, I; Z=Zvalue(cache)) for k in 1:lat.kt]
	@test relerr(g1, g2) < 1.0e-8
	g3 = cached_gf_fast(lat, K, I; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
	g3[end] = 1 - g3[1]
	@test length(g3) == lat.kτ
	# the vectorized API agrees with per-point evaluation up to implementation details
	g4 = [cached_gf(lat, (ContourIndex(k, conj=false, branch=:τ, band=1), ContourIndex(1, conj=true, branch=:τ, band=1)), K, I; cache=cache) for k in 1:lat.kτ]
	@test relerr(g3, g4) < 5.0e-3
end
