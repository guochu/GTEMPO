@testset "Few-mode bath, imaginary time（ED）" begin
	μ = 0.7; ω = 1.0; α = 0.5
	δτ = 0.1; N = 8; β = N * δτ

	# 单带单浴模（16 维 ED）
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = hybridization(lat, bath)
	model = ToulouseIM(ϵ_d=μ)
	K, Is = fermionic_setup(lat, bath, model)
	Z = integrate(K, Is...)
	# 配分函数（含浴耦合的完整 ED 谱）
	# 浴归一化（Z_B 铁律）+ 行因子化 O(δτ) 常数（~0.1/δτ，δτ=0.1）
	@test real(Z) ≈ real(tr(exp(-β * H))) / (1 + exp(-β * ω)) rtol = 4e-2
	# Matsubara G(τ, 0)（时序因子 −⟨a(τ)a†⟩）
	m = [matsubara(lat, model, kernel, Is...; i=i, j=0, Z=Z) for i in 1:lat.k]
	@test relerr(m, -g_ed[2:end]) < 2e-3
	# occupation 平衡常数
	occ = [occupation(lat, model, kernel, Is...; j=i, branch=:τ, Z=Z) for i in 1:lat.k]
	n_exact = real(tr(adag * a * exp(-β * H))) / real(tr(exp(-β * H)))
	@test maximum(abs.(occ .- n_exact)) < 1e-3
end

@testset "Few-mode bath, imaginary time: 两浴模" begin
	μ = 0.7
	δτ = 0.1; N = 10; β = N * δτ
	specs = [(1.0, 0.3), (2.0, 0.4)]
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=specs, nbands=1)
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	bath = discretefermionicbath([s[1] for s in specs], [s[2] for s in specs]; β=β)
	kernel = hybridization(lat, bath)
	model = ToulouseIM(ϵ_d=μ)
	K, Is = fermionic_setup(lat, bath, model)
	Z = integrate(K, Is...)
	m = [matsubara(lat, model, kernel, Is...; i=i, j=0, Z=Z) for i in 1:lat.k]
	@test relerr(m, -g_ed[2:end]) < 2e-3
end

@testset "Few-mode bath, imaginary time: 两带 U≠0" begin
	μ = 0.7; U = 1.0; ω = 1.0; α = 0.5
	δτ = 0.1; N = 10; β = N * δτ
	H, a, adag, H0 = singlemode_ed(μ=μ, U=U, bathspecs=[(ω, α), (ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	lat = GrassmannLattice(N=N, δτ=δτ, bands=2, contour=:imag)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = hybridization(lat, bath)
	model = AndersonIM(U=U, ϵ_d=μ)
	K, Is = fermionic_setup(lat, bath, model)
	Z = integrate(K, Is...)
	@test real(Z) ≈ real(tr(exp(-β * H))) / (1 + exp(-β * ω))^2 rtol = 4e-2
	m = [matsubara(lat, model, kernel, Is...; i=i, j=0, Z=Z) for i in 1:lat.k]
	@test relerr(m, -g_ed[2:end]) < 2e-3
	occ = occupation(lat, model, kernel, Is...; j=1, branch=:τ, band=1, Z=Z)
	n_exact = real(tr(adag * a * exp(-β * H))) / real(tr(exp(-β * H)))
	@test occ ≈ n_exact atol = 1e-3
end

@testset "Few-mode bath, imaginary time: quench（τ 支只感 hτ）" begin
	μ0 = 0.7; μ1 = -1.3; ω = 1.0; α = 0.5
	δτ = 0.1; N = 8; β = N * δτ
	H, a, adag, H0 = singlemode_ed(μ=μ0, U=0, bathspecs=[(ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	bath = discretefermionicbath([ω], [α]; β=β)
	model = QuenchedImpurityHamiltonian([tunneling(1, 1, coeff=μ0)],
	                                    [tunneling(1, 1, coeff=μ1)]; bands=1)
	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	kernel = hybridization(lat, bath)
	K, Is = fermionic_setup(lat, bath, model)
	Z = integrate(K, Is...)
	m = [matsubara(lat, model, kernel, Is...; i=i, j=0, Z=Z) for i in 1:lat.k]
	@test relerr(m, -g_ed[2:end]) < 2e-3
end
