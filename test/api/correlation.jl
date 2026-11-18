@testset "API: correlationfunction（kernel 构造，W2）" begin
	β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)

	# 虚时：τ/τ 块（CorrelationMatrix 存储）
	lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernel = correlationfunction(bath, lat)
	@test kernel isa HybridizationKernel
	@test kernel.blocks[(ContourBlockIndex(:τ, 1), ContourBlockIndex(:τ, 1))] isa CorrelationMatrix
	# 实时：四块（跨支 = CorrelationMatrix，Toeplitz 向量存储）
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	kernel = correlationfunction(bath, lat)
	@test kernel.blocks[(ContourBlockIndex(:+, 1), ContourBlockIndex(:-, 1))] isa CorrelationMatrix
	# 混合：τ↔rt = Matrix（构造期折叠物化），同支 = CorrelationMatrix
	lat = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	kernel = correlationfunction(bath, lat)
	@test kernel.blocks[(ContourBlockIndex(:τ, 1), ContourBlockIndex(:+, 1))] isa Matrix
	@test kernel.blocks[(ContourBlockIndex(:+, 1), ContourBlockIndex(:+, 1))] isa CorrelationMatrix
	# 单带构造：bands = [band]，块键仅该带
	lat2 = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed, bands=2)
	k1 = correlationfunction(bath, lat2; band=2)
	@test k1.bands == [2]
	@test keys(k1.blocks) == Set(((ContourBlockIndex(b1, 2), ContourBlockIndex(b2, 2))
	                              for b1 in branches(lat2), b2 in branches(lat2)))
	@test_throws BoundsError correlationfunction(bath, lat; band=2)   # band 越界
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
