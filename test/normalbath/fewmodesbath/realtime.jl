@testset "Few-mode bath, real time（ED）" begin
	μ = 0.7; ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 8
	ts = (1:Nt) .* δt

	# ED 参照：单浴模、可分离热初态
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = hybridization(lat, bath)
	model = ToulouseIM(ϵ_d=μ)
	for order in (1, 2)
		K, Is = fermionic_setup(lat, bath, model; β=β, order=order, coupling_norm=coupling_norm(bath))
		Z = integrate(K, Is...)
		for i in 1:Nt
			g = greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=order,
			            coupling_norm=coupling_norm(bath), initial_state=thermal_impurity_state(model, β))
			l = lesser(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=order,
			           coupling_norm=coupling_norm(bath), initial_state=thermal_impurity_state(model, β))
			@test g ≈ gt_ed[i] atol = 1e-4
			@test l ≈ lt_ed[i] atol = 1e-4
		end
	end
end

@testset "Few-mode bath, real time: 两带 U≠0" begin
	μ = 0.7; U = 1.0
	β = 1.0; δt = 0.05; Nt = 10
	ts = (1:Nt) .* δt
	ω, α = 1.0, 0.5
	H, a, adag, H0 = singlemode_ed(μ=μ, U=U, bathspecs=[(ω, α), (ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real, bands=2)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = hybridization(lat, bath; coupling_norm=coupling_norm(bath))
	model = AndersonIM(U=U, ϵ_d=μ)
	K, Is = fermionic_setup(lat, bath, model; β=β, order=2, coupling_norm=coupling_norm(bath))
	Z = integrate(K, Is...)
	for i in 1:Nt
		g = greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2,
		            coupling_norm=coupling_norm(bath), initial_state=thermal_impurity_state(model, β))
		@test g ≈ gt_ed[i] atol = 1e-4
	end
end

@testset "Few-mode bath, real time: quench Hamiltonian" begin
	μ0 = 0.7; μ1 = -0.5; ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 8
	ts = (1:Nt) .* δt
	_, a, adag, H0 = singlemode_ed(μ=μ0, U=0, bathspecs=[(ω, α)])
	H, = singlemode_ed(μ=μ1, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

	bath = discretefermionicbath([ω], [α]; β=β)
	model = QuenchedImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],
	                                    [tunneling(1, 1, coeff=μ1)]; bands=1)
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
	kernel = hybridization(lat, bath)
	Is = hybriddynamics(lat, bath, PartialIF(); band=1)
	K = sysdynamics(lat, model; initial_state=thermal_impurity_state(model, β))
	Z = integrate(K, Is)
	for i in 1:Nt
		g = greater(lat, model, kernel, Is; i=i, j=0, Z=Z, initial_state=thermal_impurity_state(model, β))
		l = lesser(lat, model, kernel, Is; i=i, j=0, Z=Z, initial_state=thermal_impurity_state(model, β))
		@test g ≈ gt_ed[i] atol = 1e-4
		@test l ≈ lt_ed[i] atol = 1e-4
	end
end

@testset "Few-mode bath, real time: 时间依赖 Hamiltonian" begin
	μ0 = 0.7; μ1 = -0.5; A = 0.5; ωt = 1.0
	ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 8

	# ED 参照：逐步常数 H 离散化（与 GTEMPO 同步长），只有截断误差
	_, a, adag, H0 = singlemode_ed(μ=μ0, U=0, bathspecs=[(ω, α)])
	ρ = exp(-β * H0); ρ /= tr(ρ)
	gt_ed = complex.(zeros(Nt))
	lt_ed = complex.(zeros(Nt))
	U = one(H0)
	for k in 1:Nt
		Hk, = singlemode_ed(μ=μ1 + A * sin(ωt * (k - 1) * δt), U=0, bathspecs=[(ω, α)])
		U = exp(-im * δt * Hk) * U
		gt_ed[k] = -im * tr((U' * a * U) * adag * ρ) / tr(ρ)
		lt_ed[k] = im * tr(adag * (U' * a * U) * ρ) / tr(ρ)
	end

	bath = discretefermionicbath([ω], [α]; β=β)
	model = TdImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],
	                              [tunneling(1, 1, coeff=μ1)],
	                              [TdImpurityOp([tunneling(1, 1)], t -> A * sin(ωt * t))])
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
	kernel = hybridization(lat, bath)
	Is = hybriddynamics(lat, bath, PartialIF(); band=1)
	# 热初态取 hτ（预淬 μ0）
	ρ0 = Matrix{ComplexF64}(fock_thermalstate(QuenchedImpurityHamiltonian(
		[tunneling(1, 1, coeff=μ0)], [tunneling(1, 1, coeff=μ1)]; bands=1), β).data)
	K = sysdynamics(lat, model; initial_state=ρ0)
	Z = integrate(K, Is)
	for i in 1:Nt
		g = greater(lat, model, kernel, Is; i=i, j=0, Z=Z, initial_state=ρ0)
		# Td 逐步常数 H 合成的 O(δt²) 相位（两侧同离散化、阶累积差异）
		@test g ≈ gt_ed[i] atol = 1e-3
	end
end
