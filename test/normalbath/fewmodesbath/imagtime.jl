@testset "Few-mode bath, imaginary time: orderings & algorithms" begin
	μ = 0.7; ω = 1.0; α = 0.5
	δτ = 0.1; N = 8; β = N * δτ
	trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

	# ED reference: one fermionic bath mode
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	spec = DiracDelta(ω=ω, α=α)
	bath = fermionicbath(spec, β=β)
	# analytic Toulouse solution for the δ spectrum — validates the ED conventions
	g_analytic = [toulouse_Gτ(bath, τ; ϵ_d=μ) for τ in 0:δτ:β]
	@test relerr(g_ed, g_analytic) < 1.0e-2

	model = ToulouseIM(μ=μ)
	for ordering in imag_orderings
		lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag, ordering=ordering)
		corr = correlationfunction(bath, lat)
		# all algorithms on the default ordering, the default algorithm on all orderings
		algs = (ordering == imag_orderings[1]) ? if_algs(trunc) :
			   [("PartialIF", PartialIF(trunc=trunc))]
		for (name, alg) in algs
			mpsI = hybriddynamics(lat, corr, alg)
			mpsK = sysdynamics(lat, model, trunc=trunc)
			mpsK = boundarycondition!(mpsK, lat)
			g = gτ_series(lat, mpsK, mpsI)
			@test relerr(g, g_ed) < 1.0e-2
		end
	end
end

@testset "Few-mode bath, imaginary time: two bath modes" begin
	μ = 0.7
	δτ = 0.1; N = 10; β = N * δτ
	trunc = truncdimcutoff(D=150, ϵ=1.0e-10)
	specs = [(1.0, 0.3), (2.0, 0.4)]

	# single-band impurity coupled to two bath modes (8-dimensional ED)
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=specs, nbands=1)
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	# the two δ modes are combined into one DiscreteSpectrum bath:
	# DiscreteSpectrum(ws, fs) represents δ peaks at ws with weights fs
	spec = DiscreteSpectrum([s[1] for s in specs], [s[2] for s in specs])
	corr = correlationfunction(fermionicbath(spec, β=β), lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, ToulouseIM(μ=μ), trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	g = gτ_series(lat, mpsK, mpsI)
	@test relerr(g, g_ed) < 2.0e-2
end

@testset "Few-mode bath, imaginary time: two bands with U" begin
	μ = 0.7; U = 1.0; ω = 1.0; α = 0.5
	δτ = 0.1; N = 10; β = N * δτ
	trunc = truncdimcutoff(D=80, ϵ=1.0e-10)

	# one bath mode per band
	H, a, adag, H0 = singlemode_ed(μ=μ, U=U, bathspecs=[(ω, α), (ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag, bands=2)
	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = AndersonIM(U=U, μ=μ)
	mpsK, Is = fermionic_setup(lat, bath, model, trunc)
	cache = environments(lat, mpsK, Is...)
	g = cached_gf_fast(lat, mpsK, Is...; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
	@test relerr(g, g_ed) < 1.0e-2

	# occupation number ⟨n₁⟩ as an additional observable
	n_gtempo = real.([occupation(lat, i, mpsK, Is..., Z=Zvalue(cache), band=1) for i in 1:lat.k])
	ρ = exp(-β * H) / tr(exp(-β * H))
	n_exact = real(tr(ρ * adag * a))
	@test abs(n_gtempo[1] - n_exact) < 1.0e-2
end

@testset "Few-mode bath, imaginary time: quench Hamiltonian" begin
	# on the imaginary axis only the pre-quench h0 enters (the :τ branch
	# propagator); h1 is deliberately very different and must not matter
	μ0 = 0.7; μ1 = -1.3; ω = 1.0; α = 0.5
	δτ = 0.1; N = 8; β = N * δτ
	trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

	H, a, adag, H0 = singlemode_ed(μ=μ0, U=0, bathspecs=[(ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = QuenchedImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],
										[tunneling(1, 1, coeff=μ1)]; bands=1)
	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	corr = correlationfunction(bath, lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	g = gτ_series(lat, mpsK, mpsI)
	@test relerr(g, g_ed) < 1.0e-2
end

@testset "Few-mode bath, imaginary time: in-place ExactTTIIF hybriddynamics!" begin
	# the whole workflow on one state: sysdynamics builds K, then
	# hybriddynamics! merges the influence functional into K in place
	μ = 0.7; ω = 1.0; α = 0.5
	δτ = 0.1; N = 8; β = N * δτ
	trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = ToulouseIM(μ=μ)
	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	corr = correlationfunction(bath, lat)
	alg = ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0)

	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	mpsK = hybriddynamics!(mpsK, lat, corr, alg)
	g = gτ_series(lat, mpsK)
	@test relerr(g, g_ed) < 1.0e-2

	# agrees with the standard separate-construction workflow
	mpsI = hybriddynamics(lat, corr, alg)
	mpsK0 = sysdynamics(lat, model, trunc=trunc)
	mpsK0 = boundarycondition!(mpsK0, lat)
	g0 = gτ_series(lat, mpsK0, mpsI)
	@test relerr(g, g0) < 1.0e-3
end

@testset "Few-mode bath, imaginary time: time-dependent Hamiltonian" begin
	# on the imaginary axis only the constant hτ enters; ht/htt are
	# deliberately very different and must not matter
	μ0 = 0.7; μ1 = -1.3; A = 5.0; ωt = 1.0
	ω = 1.0; α = 0.5
	δτ = 0.1; N = 8; β = N * δτ
	trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

	H, a, adag, H0 = singlemode_ed(μ=μ0, U=0, bathspecs=[(ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = TdImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],
									[tunneling(1, 1, coeff=μ1)],
									[TdImpurityOp([tunneling(1, 1)], t -> A * sin(ωt * t))])
	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	corr = correlationfunction(bath, lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	g = gτ_series(lat, mpsK, mpsI)
	@test relerr(g, g_ed) < 1.0e-2
end
