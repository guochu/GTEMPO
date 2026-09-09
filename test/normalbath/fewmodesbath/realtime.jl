@testset "Few-mode bath, real time: orderings & algorithms" begin
	μ = 0.7; ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 8
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

	# ED reference: one fermionic bath mode, separable thermal initial state
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

	spec = DiracDelta(ω=ω, α=α)
	bath = fermionicbath(spec, β=β)
	model = ToulouseIM(μ=μ)
	for ordering in real_orderings
		lat = GrassmannLattice(N=Nt, δt=δt, contour=:real, ordering=ordering)
		corr = correlationfunction(bath, lat)
		# all algorithms on the default ordering, the default algorithm on all orderings
		algs = (ordering == real_orderings[1]) ? if_algs(trunc) :
			   [("PartialIF", PartialIF(trunc=trunc))]
		for (name, alg) in algs
			mpsI = hybriddynamics(lat, corr, alg)
			mpsK = sysdynamics(lat, model, trunc=trunc)
			mpsK = boundarycondition!(mpsK, lat)
			mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)
			gt, lt = gtlt_series(lat, mpsK, mpsI)
			@test relerr(gt, gt_ed) < 3.0e-2
			@test relerr(lt, lt_ed) < 3.0e-2
		end
	end
end

@testset "Few-mode bath, real time: two bands with U" begin
	μ = 0.7; U = 1.0
	β = 1.0; δt = 0.05; Nt = 10
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=80, ϵ=1.0e-10)
	ω, α = 1.0, 0.5

	# two bands with U: ED reference (16-dimensional ED), one bath mode per band
	H, a, adag, H0 = singlemode_ed(μ=μ, U=U, bathspecs=[(ω, α), (ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real, bands=2)
	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	mpsK, Is = fermionic_setup(lat, bath, AndersonIM(U=U, μ=μ), trunc, β=β)
	gt, lt = gtlt_series(lat, mpsK, Is...)
	@test relerr(gt, gt_ed) < 3.0e-2
	@test relerr(lt, lt_ed) < 3.0e-2
end

@testset "Few-mode bath, real time: quench Hamiltonian" begin
	# the impurity is thermalized with the pre-quench h0 (μ0, also used by the
	# separable thermal initial state) and evolves with h1 (μ1) for t > 0
	μ0 = 0.7; μ1 = -0.5; ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 8
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

	# ED: separable initial state from the decoupled μ0 Hamiltonian, real-time
	# evolution under the coupled μ1 Hamiltonian
	_, a, adag, H0 = singlemode_ed(μ=μ0, U=0, bathspecs=[(ω, α)])
	H, = singlemode_ed(μ=μ1, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = QuenchedImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],
										[tunneling(1, 1, coeff=μ1)]; bands=1)
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
	corr = correlationfunction(bath, lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)
	gt, lt = gtlt_series(lat, mpsK, mpsI)
	@test relerr(gt, gt_ed) < 3.0e-2
	@test relerr(lt, lt_ed) < 3.0e-2
end

@testset "Few-mode bath, real time: time-dependent Hamiltonian" begin
	# the impurity level is modulated as μ(t) = μ1 + A·sin(ωt·t): the thermal
	# initial state is built from hτ (μ0), the real-time branches evolve with
	# H(t). The ED reference uses the same per-step constant-H discretization
	# as GTEMPO, so only truncation errors remain
	μ0 = 0.7; μ1 = -0.5; A = 0.5; ωt = 1.0
	ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 8
	trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

	_, a, adag, H0 = singlemode_ed(μ=μ0, U=0, bathspecs=[(ω, α)])
	ρ = exp(-β * H0); ρ /= tr(ρ)
	gt_ed = complex.(zeros(Nt+1))
	lt_ed = complex.(zeros(Nt+1))
	U = one(H0)
	for k in 0:Nt
		if k > 0
			Hk, = singlemode_ed(μ=μ1 + A * sin(ωt * (k-1) * δt), U=0, bathspecs=[(ω, α)])
			U = exp(-im * δt * Hk) * U
		end
		# U accumulates the forward evolution e^{-iδtH}; the Heisenberg
		# operator is a(t) = U' * a * U. G^>(t, 0) = -i⟨a(t) a†⟩,
		# G^<(t, 0) = i⟨a† a(t)⟩ (the creation operator sits at t = 0)
		gt_ed[k+1] = -im * tr((U' * a * U) * adag * ρ)
		lt_ed[k+1] = im * tr(adag * (U' * a * U) * ρ)
	end

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = TdImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],
									[tunneling(1, 1, coeff=μ1)],
									[TdImpurityOp([tunneling(1, 1)], t -> A * sin(ωt * t))])
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
	corr = correlationfunction(bath, lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)
	gt, lt = gtlt_series(lat, mpsK, mpsI)
	@test relerr(gt, gt_ed) < 3.0e-2
	@test relerr(lt, lt_ed) < 3.0e-2
end
