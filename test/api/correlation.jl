@testset "API: correlationfunction（kernel 构造，W2）" begin
	β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)

	# 虚时：τ/τ 块（CorrelationMatrix 存储）
	lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernel = correlationfunction(bath, lat)
	@test kernel isa HybridizationKernel
	@test kernel.blocks[(ContourBlockIndex(:τ, 1, true), ContourBlockIndex(:τ, 1, false))] isa CorrelationMatrix
	@test eltype(first(values(kernel.blocks))) == Float64   # 实谱虚时核 = Float64（T = scalartype(corr)）
	# 实时：四块（跨支 = CorrelationMatrix，Toeplitz 向量存储）
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	kernel = correlationfunction(bath, lat)
	@test kernel.blocks[(ContourBlockIndex(:+, 1, true), ContourBlockIndex(:-, 1, false))] isa CorrelationMatrix
	# 混合：τ↔rt = Matrix（构造期折叠物化），同支 = CorrelationMatrix
	lat = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	kernel = correlationfunction(bath, lat)
	@test kernel.blocks[(ContourBlockIndex(:τ, 1, true), ContourBlockIndex(:+, 1, false))] isa Matrix
	@test kernel.blocks[(ContourBlockIndex(:+, 1, true), ContourBlockIndex(:+, 1, false))] isa CorrelationMatrix
	# 单带构造：bands = [band]，块键仅该带
	lat2 = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed, bands=2)
	k1 = correlationfunction(bath, lat2; band=2)
	@test k1.bands == [2]
	@test keys(k1.blocks) == Set(((ContourBlockIndex(b1, 2, true), ContourBlockIndex(b2, 2, false))
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

@testset "API: kernel 动态标量型（实谱 Float64 / 跨标量相加）" begin
	β = 1.0
	# 实谱虚时核 = Float64（T = scalartype(corr)，lattice scalartype Float64 同型）
	bath = discretefermionicbath([1.0], [0.5]; β=β)
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kimg = correlationfunction(bath, latτ)
	@test eltype(first(values(kimg.blocks))) == Float64
	# 跨标量相加：同键 = QuAPI CM 向量和（promote 复）；缺键 = 透传标量对齐 T
	blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{ComplexF64},CorrelationMatrix{ComplexF64}}}()
	key = (ContourBlockIndex(:τ, 1, true), ContourBlockIndex(:τ, 1, false))
	blocks[key] = CorrelationMatrix(ComplexF64[0.3, 0.2, 0.1, 0.05], ComplexF64[0.3, 0.2, 0.1, 0.05])
	kcx = HybridizationKernel(blocks, 0.4, [1])
	ksum = kimg + kcx
	@test ksum.blocks[key] isa CorrelationMatrix
	@test eltype(ksum.blocks[key]) == ComplexF64
	@test Matrix(ksum.blocks[key]) ≈ Matrix(kimg.blocks[key]) + Matrix(kcx.blocks[key]) atol = 1e-14
	@test ksum.coupling_norm ≈ kimg.coupling_norm + 0.4
	# 核相加的 CM 保持（同标量 CM+CM 不物化稠密）
	ka = correlationfunction(discretefermionicbath([2.0], [0.3]; β=β), latτ)
	kab = kimg + ka
	@test kab.blocks[key] isa CorrelationMatrix
	@test Matrix(kab.blocks[key]) ≈ Matrix(kimg.blocks[key]) + Matrix(ka.blocks[key]) atol = 0.0
end

@testset "API: scheme 面（IFD 后端切换 + 严格 IF 入口）" begin
	import IFDiscretization as IFD
	# plan/ifd_backend_switch：后端切换 = correlationfunction 第三位 scheme 参数
	# （缺省 IFD.QuAPI()）；严格面 = imag + Exact 元组解包 + coupling_norm 复用（= normdelta）。
	β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	latr = GrassmannLattice(N=4, δt=0.1, contour=:real)
	latm = GrassmannLattice(Nt=2, δt=0.1, Nτ=3, δτ=0.25, contour=:mixed)

	# T3a：缺省 ≡ 显式 QuAPI scheme（逐位中性证明）
	k1 = correlationfunction(bath, latτ)
	k2 = correlationfunction(bath, latτ, IFD.QuAPI())
	@test collect(pairs(k1.blocks)) == collect(pairs(k2.blocks))


	# T3b：imag + Exact（求解器要求可调用谱 → 连续谱 bath）——blocks ≡ 元组
	# 第一位逐位、normdelta ≡ 第二位（det(BB)）
	spec = spectrum_func()
	bathc = fermionicbath(spec, β=β, μ=0.0)
	scheme = IFD.ExactEqualDistantKernel(IFD.EDToulouseSolver(8))
	kx = correlationfunction(bathc, latτ, scheme)
	res = IFD.Δτ(bathc; N=4, δτ=0.25, scheme)
	@test kx.coupling_norm == real(res[2])   # coupling_norm 语义随 scheme（exact = normdelta 透传）；虚部尾守卫取实
	@test Matrix(kx.blocks[(ContourBlockIndex(:τ, 1, true), ContourBlockIndex(:τ, 1, false))]) ==
	      Matrix(res[1].data[1:4, 1:4])

	# T4：real/mixed + Exact 拒绝（strict 容器 = 装配级骨架+踢契约，非朴素指数
	# 系数——消费路线归 exactif_port_plan）；非法 scheme → MethodError 自然冒泡
	@test_throws ArgumentError correlationfunction(bathc, latr, scheme)
	@test_throws ArgumentError correlationfunction(bathc, latm, scheme)
	@test_throws ArgumentError correlationfunction(bath, latr, scheme)
	@test_throws MethodError correlationfunction(bath, latτ, "bogus")

	# T5：exact IF 装配 smoke（normdelta·Pf(B) 的 Z 消费契约归 exactif_port_plan
	# ——本计划只保证字段忠实透传 + 可装配：exact IF 的 bath 归一 Z 与 QuAPI
	# scheme 同量级）
	Isq, _ = hybriddynamics(latτ, k1, PartialIF())
	Ise, _ = hybriddynamics(latτ, kx, PartialIF())
	model = ToulouseIM(ϵ_d=1.25 * pi)
	K = sysdynamics(latτ, model)
	Zq = integrate(K, Isq)
	Zx = integrate(K, Ise)
	@test isfinite(Zx) && 0.5 < abs(Zx) / abs(Zq) < 2.0
end
