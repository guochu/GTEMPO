@testset "Continuous bath, mixed time（Kadanoff）" begin
	β = 1.0; δt = 0.05; Nt = 3; δτ = 0.25; Nτ = 4
	bath = fermionicbath(spectrum_func(), β=β, μ=0)
	lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)
	kernel = correlationfunction(bath, lat)
	Is, _ = hybriddynamics(lat, kernel, PartialIF())
	model = ToulouseIM(ϵ_d=1.0)
	K = sysdynamics(lat, model)
	Z = integrate(K, Is)
	# 结构性：各观测量有限、τ 支 occupation 平衡恒定
	occs = [occupation(lat, model, kernel, Is; j=i, branch=:τ, Z=Z) for i in 1:Nτ]
	# Kadanoff τ 腿离散化 → 占据的 O(δτ) 平衡漂移
	@test maximum(abs.(occs .- occs[1])) < 5e-3
	for i in 1:Nt
		g = greater(lat, model, kernel, Is; i=i, j=0, Z=Z)
		@test isfinite(real(g)) && isfinite(imag(g))
		l = lesser(lat, model, kernel, Is; i=i, j=0, Z=Z)
		@test isfinite(real(l)) && isfinite(imag(l))
	end
	# τ 支 matsubara 与纯虚时轮廓一致（同一 β 谱）
	latτ = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag)
	kernelτ = correlationfunction(bath, latτ)
	Isτ, _ = hybriddynamics(latτ, kernelτ, PartialIF())
	Kτ = sysdynamics(latτ, model)
	Zτ = integrate(Kτ, Isτ)
	for i in 1:Nτ
		m1 = matsubara(lat, model, kernel, Is; i=i, j=0, Z=Z)
		m2 = matsubara(latτ, model, kernelτ, Isτ; i=i, j=0, Z=Zτ)
		# 混合 vs 纯虚时 τ 支：两构造的一阶 IF 常数不同 → O(δτ) 差
		@test m1 ≈ m2 atol = 5e-3
	end
end
