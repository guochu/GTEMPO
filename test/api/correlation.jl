@testset "API: correlationfunction（QuAPI 包装）" begin
	β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)

	# 虚时：Δτ 矩阵
	lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	corr = correlationfunction(bath, lat)
	@test corr isa ImagCorrelationFunction
	# 实时：QuAPI 关联对象四块
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	corr = correlationfunction(bath, lat)
	@test corr isa RealCorrelationFunction
	# 混合
	lat = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	corr = correlationfunction(bath, lat)
	@test corr isa MixedCorrelationFunction
end

@testset "API: baths（ImpurityModelBase 皮肤）" begin
	β = 1.0
	bath = discretefermionicbath([1.0, 2.0], [0.5, 0.3]; β=β)
	@test bath.β == β
	bath = fermionicbath(spectrum_func(), β=β, μ=0.3)
	@test bath.β == β && bath.μ == 0.3
	@test bath isa AbstractNormalBath
end

@testset "API: coupling_norm（浴谱总强度）" begin
	β = 1.0
	# δ 谱：α；离散谱：Σαₙ（多点）
	@test coupling_norm(DiracDelta(ω=1.0, α=0.5)) == 0.5
	ds = DiscreteSpectrum([1.0, 2.0, 3.0], [0.5, 0.3, 0.2])
	@test coupling_norm(ds) == 1.0
	# 连续谱：矩形解析对照 + 高斯谱（8σ 截断，尾质量 < 1e-14）
	spec_rect = spectrum(ϵ -> 0.5, lb=-1.0, ub=2.0)
	@test coupling_norm(spec_rect) ≈ 0.5 * 3.0 rtol = 1e-10
	spec_g = spectrum(ϵ -> exp(-(ϵ - 0.3)^2 / 0.02), lb=0.3 - 0.8, ub=0.3 + 0.8)
	@test coupling_norm(spec_g) ≈ sqrt(0.02 * π) * 0.9999999999 rtol = 1e-6
	# δ×函数乘积（DeltaMultF/DeltasMultF；g ≡ 1 → 与直接公式恒等；加权抽查）
	@test coupling_norm(DiracDelta(ω=1.0, α=0.5) * (ϵ -> 2.0)) ≈ 1.0
	@test coupling_norm(ds * (ϵ -> 2.0)) ≈ 2.0
	@test coupling_norm(DiracDelta(ω=1.0, α=0.5) * (ϵ -> ϵ)) ≈ 0.5
	# bath 层：谱直取——μ/β 无关
	spec = spectrum_func()
	@test coupling_norm(fermionicbath(spec, β=β, μ=0.0)) ≈ coupling_norm(fermionicbath(spec, β=2.0, μ=0.3))
	@test coupling_norm(discretefermionicbath([1.0, 2.0], [0.5, 0.3]; β=β)) == 0.8
	@test coupling_norm(fermionicbath(spec, β=β, μ=0.0)) ≈ 5.0 atol = 1e-6      # 半圆谱 [-10,10]：∫ sqrt(1-(ε/D)²)/π dε = 10·(π/2)/π = 5
end
