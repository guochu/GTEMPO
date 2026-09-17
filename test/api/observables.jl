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

@testset "API: greater / lesser / contour ordered gf" begin
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	model = ToulouseIM(μ=-0.5)
	bath = fermionicbath(spectrum_func(), β=1.0, μ=0)
	lat = GrassmannLattice(N=4, δt=0.05, contour=:real)
	corr = correlationfunction(bath, lat)
	I = hybriddynamics(lat, corr, trunc=trunc)
	K = boundarycondition!(sysdynamics(lat, model, trunc=trunc), lat)
	Z = integrate(lat, K, I)
	cache = environments(lat, K, I)

	# greater ⟨aᵢ bⱼ⟩ (i ≥ j, both variables on the + branch)
	g1 = greater(lat, 3, 1, K, I; Z=Z)
	@test g1 == gf(lat, (ContourIndex(3, conj=false, branch=:+, band=1), ContourIndex(1, conj=true, branch=:+, band=1)), K, I; Z=Z)
	# single-argument version defaults to j = 1
	@test greater(lat, 3, K, I; Z=Z) == g1
	# lesser requires i ≤ j
	l1 = lesser(lat, 1, 3, K, I; Z=Z)
	@test l1 == gf(lat, (ContourIndex(1, conj=true, branch=:-, band=1), ContourIndex(3, conj=false, branch=:+, band=1)), K, I; Z=Z)
	@test lesser(lat, 3, K, I; Z=Z) == l1
	# contour ordering: a < b gives -(b, a), otherwise (a, b); conj(a)=false, conj(b)=true
	# (contour_ordered_gf requires a real Z, unlike greater/lesser)
	x = ContourIndex(1, conj=false, branch=:+, band=1)
	y = ContourIndex(3, conj=true, branch=:+, band=1)
	@test contour_ordered_gf(lat, x, y, K, I; Z=real(Z)) ≈ -gf(lat, (y, x), K, I; Z=Z) atol=1.0e-9
	x2 = ContourIndex(3, conj=false, branch=:+, band=1)
	y2 = ContourIndex(1, conj=true, branch=:+, band=1)
	@test contour_ordered_gf(lat, x2, y2, K, I; Z=real(Z)) ≈ g1 atol=1.0e-9
	@test_throws ArgumentError contour_ordered_gf(lat, y2, x2, K, I; Z=1.0)
	# cached versions agree with the direct evaluation (up to roundoff)
	@test cached_greater(lat, 3, 1, K, I; cache=cache) ≈ g1 atol=1.0e-12
	@test cached_lesser(lat, 1, 3, K, I; cache=cache) ≈ l1 atol=1.0e-12
	cc1 = cached_contour_ordered_gf(lat, x, y, K, I; cache=cache)
	cc2 = cached_contour_ordered_gf(lat, x2, y2, K, I; cache=cache)
	# the a > b branch returns the greater correlation directly, the a < b
	# branch returns -⟨b a⟩ (contour-ordered sign)
	@test cc2 ≈ cached_greater(lat, 3, 1, K, I; cache=cache) atol=1.0e-12
	@test cc1 ≈ -cached_gf(lat, (y, x), K, I; cache=cache) atol=1.0e-12
	# vectorized fast variants
	cg = cached_greater_fast(lat, K, I; cache=cache)
	cl = cached_lesser_fast(lat, K, I; cache=cache)
	@test length(cg) == length(cl) == lat.k
	@test cg[3] ≈ g1 atol=1.0e-8
	@test cl[2] ≈ cached_lesser(lat, 1, 2, K, I; cache=cache) atol=1.0e-8
end

@testset "API: electric current & heat current" begin
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	model = ToulouseIM(μ=-0.5)
	bath = fermionicbath(spectrum_func(), β=1.0, μ=0)
	lat = GrassmannLattice(N=4, δt=0.05, contour=:real)
	corr = correlationfunction(bath, lat)
	I = hybriddynamics(lat, corr, trunc=trunc)
	K = boundarycondition!(sysdynamics(lat, model, trunc=trunc), lat)
	Z = integrate(lat, K, I)
	cache = environments(lat, K, I)

	# electric current: vector over k = 2:lattice.k, scalar at a fixed k
	J = electriccurrent(lat, corr, K, I; Z=Z)
	@test J isa Vector && length(J) == lat.k - 1
	@test electriccurrent(lat, corr, lat.k, K, I; Z=Z) ≈ J[end] atol=1.0e-10
	# MPO-based fast variant
	Jf = electriccurrent_fast(lat, corr, K, I)
	@test length(Jf) == length(J) && relerr(Jf, J) < 1.0e-6
	@test electriccurrent_fast(lat, corr, lat.k, K, I) ≈ Jf[end] atol=1.0e-10
	# cached variants
	Jc = cached_electriccurrent(lat, corr, K, I; cache=cache)
	Jcf = cached_electriccurrent_fast(lat, corr, K, I; cache=cache)
	@test relerr(Jc, J) < 1.0e-8
	@test relerr(Jcf, Jf) < 1.0e-8

	# heat current: the bath spectrum weighted by ω, otherwise identical to the current
	hcorr = heatcorrelationfunction(bath, lat)
	@test hcorr isa RealCorrelationFunction && size(hcorr.G₊₊) == size(corr.G₊₊)
	Jh = heatcurrent_fast(lat, bath, K, I)
	@test relerr(Jh, electriccurrent_fast(lat, hcorr, K, I)) == 0
	Jhc = cached_heatcurrent_fast(lat, bath, K, I; cache=cache)
	@test relerr(Jhc, cached_electriccurrent_fast(lat, hcorr, K, I; cache=cache)) == 0
end

@testset "API: nn family & expectationvalue" begin
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	model = ToulouseIM(μ=-0.5)
	bath = fermionicbath(spectrum_func(), β=1.0, μ=0)
	lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	corr = correlationfunction(bath, lat)
	I = hybriddynamics(lat, corr, trunc=trunc)
	K = boundarycondition!(sysdynamics(lat, model, trunc=trunc), lat)
	Z = integrate(lat, K, I)
	cache = environments(lat, K, I)

	# ⟨n̂ᵢ n̂ⱼ⟩: two independent implementations (operator insertion and the
	# four-point Green's function), each checked for index symmetry
	v = nn(lat, 2, 3, K, I; Z=Z)
	v2 = nn2(lat, 2, 3, K, I; Z=Z)
	@test isreal(v) && isreal(v2)
	@test v ≈ nn(lat, 3, 2, K, I; Z=Z) atol=1.0e-10    # n̂ᵢ n̂ⱼ = n̂ⱼ n̂ᵢ
	@test v2 ≈ nn2(lat, 3, 2, K, I; Z=Z) atol=1.0e-10
	# diagonal element equals the occupation
	@test nn(lat, 2, 2, K, I; Z=Z) ≈ occupation(lat, 2, K, I; Z=Z) atol=1.0e-8
	# insert_n: copy and in-place versions agree
	K2 = insert_n(lat, K, 2)
	@test distance(K2, insert_n!(lat, deepcopy(K), 2)) == 0
	# the n̂-inserted state reproduces ⟨n̂₂⟩
	@test relerr(integrate(lat, K2, I) / Z, nn(lat, 2, 2, K, I; Z=Z)) < 1.0e-8
	# cached versions
	@test cached_nn(lat, 2, 3, K, I; cache=cache) ≈ v atol=1.0e-8
	@test cached_nn2(lat, 2, 3, K, I; cache=cache) ≈ v2 atol=1.0e-8
	# expectationvalue of the two-point GTerm ⟨d₂ d†₁⟩ from the cache
	a = ContourIndex(2, conj=false, branch=:τ, band=1)
	b = ContourIndex(1, conj=true, branch=:τ, band=1)
	gref = gf(lat, (a, b), K, I; Z=Z)
	@test expectationvalue(GTerm(lat[a], lat[b], coeff=1), cache) ≈ gref atol=1.0e-10
	@test expectationvalue(convert(PartialMPO, GTerm(lat[a], lat[b], coeff=1)), cache) ≈ gref atol=1.0e-10
end
