@testset "Few-mode bath, mixed time: orderings" begin
	μ = 0.7; ω = 1.0; α = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=80, ϵ=1.0e-10)

	# ED references with the interacting thermal state (built by the τ leg)
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
	gτ_ref = gτ_ed(H, a, adag, 0:δτ:β, β)

	spec = DiracDelta(ω=ω, α=α)
	bath = fermionicbath(spec, β=β)
	model = ToulouseIM(μ=μ)
	for ordering in mixed_orderings
		lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, ordering=ordering)
		corr = correlationfunction(bath, lat)
		mpsI = hybriddynamics(lat, corr, trunc=trunc)
		mpsK = sysdynamics(lat, model, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lat)
		gt, lt, gτ = gtltgτ_series(lat, mpsK, mpsI)
		@test relerr(gt, gt_ed) < 3.0e-2
		@test relerr(lt, lt_ed) < 3.0e-2
		@test relerr(gτ, gτ_ref) < 1.0e-2
	end
end

@testset "Few-mode bath, mixed time: two bands with U" begin
	μ = 0.7; U = 1.0; ω = 1.0; α = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)

	# one bath mode per band
	H, a, adag, H0 = singlemode_ed(μ=μ, U=U, bathspecs=[(ω, α), (ω, α)])
	gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
	gτ_ref = gτ_ed(H, a, adag, 0:δτ:β, β)

	lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, bands=2)
	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	mpsK, Is = fermionic_setup(lat, bath, AndersonIM(U=U, μ=μ), trunc)
	gt, lt, gτ = gtltgτ_series(lat, mpsK, Is...)
	@test relerr(gt, gt_ed) < 3.0e-2
	@test relerr(lt, lt_ed) < 3.0e-2
	@test relerr(gτ, gτ_ref) < 1.0e-2
end

@testset "Few-mode bath, mixed time: two-band interacting mixed GF (full grid)" begin
	mu = 0.7; U = 1.0; w = 1.0; alpha = 0.5
	beta = 1.0; dtau = 0.1; Ntau = 10
	dt = 0.05; Nt = 6
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)

	# ED Hamiltonian: 2 impurity bands + 2 bath modes (dim = 16, exact reference
	# including the U interaction). The mixed GF G(τ_i, t_j) = <d(τ_i)d†(t_j)> is
	# checked on the full (τ, t) grid for both real branches.
	H, a, adag, H0 = singlemode_ed(μ=mu, U=U, bathspecs=[(w, alpha), (w, alpha)])
	rho = exp(-beta * H); rho /= tr(rho)

	spec = DiracDelta(ω=w, α=alpha)
	bath = fermionicbath(spec, β=beta)
	model = AndersonIM(U=U, μ=mu)

	lat = GrassmannLattice(Nt=Nt, δt=dt, Nτ=Ntau, δτ=dtau, contour=:mixed, bands=2)
	mpsK, Is = fermionic_setup(lat, bath, model, trunc)
	cache = environments(lat, mpsK, Is...)

	function mixed_gtempo(i, br, j)
		cached_gf(lat, (ContourIndex(i, conj=false, branch=:τ, band=1),
		                ContourIndex(j, conj=true, branch=br, band=1)), mpsK, Is...; cache=cache)
	end
	function mixed_ed(tau_v, t_v)
		op1 = exp(tau_v * H) * a * exp(-tau_v * H)
		op2 = exp(im * t_v * H) * adag * exp(-im * t_v * H)
		tr(rho * op1 * op2) / tr(rho)
	end

	tol = 3.0e-2
	for br in (:+, :-)
		for i in 1:Ntau+1
			tau_v = (i - 1) * dtau
			for j in 1:Nt+1
				@test _error(mixed_gtempo(i, br, j), mixed_ed(tau_v, (j - 1) * dt), tol) < 1
			end
		end
	end
end

@testset "Few-mode bath, mixed time: quench Hamiltonian" begin
	# Kadanoff-contour quench: the τ leg (with the bath hybridization on the
	# τ branch) builds the interacting thermal state e^{-βH0} of the coupled
	# pre-quench system (h0, μ0), the real-time branches evolve with h1 (μ1)
	μ0 = 0.7; μ1 = -0.5; ω = 1.0; α = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=80, ϵ=1.0e-10)

	# ED: initial state exp(-β H0c) of the coupled μ0 Hamiltonian, real-time
	# evolution under the coupled μ1 Hamiltonian
	H0c, a, adag, _ = singlemode_ed(μ=μ0, U=0, bathspecs=[(ω, α)])
	H1c, = singlemode_ed(μ=μ1, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H1c, a, adag, H0c, ts, β)
	gτ_ref = gτ_ed(H0c, a, adag, 0:δτ:β, β)

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = QuenchImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],
										[tunneling(1, 1, coeff=μ1)]; bands=1)
	lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)
	corr = correlationfunction(bath, lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	gt, lt, gτ = gtltgτ_series(lat, mpsK, mpsI)
	@test relerr(gt, gt_ed) < 3.0e-2
	@test relerr(lt, lt_ed) < 3.0e-2
	@test relerr(gτ, gτ_ref) < 1.0e-2
end
