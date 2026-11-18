@testset "Few-mode bath, mixed time（Kadanoff; ED）" begin
	μ = 0.7; ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 4; δτ = 0.25; Nτ = 4
	ts = (1:Nt) .* δt

	# ED 参照：τ 支构建相互作用热态 exp(−βH)
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)

	lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = correlationfunction(bath, lat)
	model = ToulouseIM(ϵ_d=μ)
	Is, _ = hybriddynamics(lat, kernel, PartialIF())
	K = sysdynamics(lat, model)
	Z = integrate(K, Is)
	for i in 1:Nt
		g = greater(lat, model, kernel, Is; i=i, j=0, Z=Z)
		l = lesser(lat, model, kernel, Is; i=i, j=0, Z=Z)
		@test g ≈ gt_ed[i] atol = 5e-3
		@test l ≈ lt_ed[i] atol = 5e-3
	end
	# 虚时支：matsubara 同 ED 热态
	τs = (1:Nτ) .* δτ
	gτ_ref = correlation_2op_1τ(H, a, adag, τs; β=β)
	for i in 1:Nτ
		m = matsubara(lat, model, kernel, Is; i=i, j=0, Z=Z)
		@test m ≈ -gτ_ref[i] atol = 2e-3
	end
end
