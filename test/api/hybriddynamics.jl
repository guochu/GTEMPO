using LinearAlgebra: norm
β0 = 1.0
bath0 = discretefermionicbath([1.0], [0.5]; β=β0)

@testset "API: hybridization kernel" begin
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	kernel = hybridization(lat, bath0; coupling_norm=0.5)
	@test kernel isa HybridizationKernel
	for b1 in (:+, :-), b2 in (:+, :-)
		@test size(kernel.blocks[(b1, b2)]) == (4, 4)
	end
	@test kernel.coupling_norm == 0.5
	@test_throws ArgumentError hybridization(lat, bath0; coupling_norm=-1.0)
	# 形状守卫（构造端预建校验）：corr 窗口尺寸须覆盖 lattice 网格
	lat2 = GrassmannLattice(N=6, δt=0.1, contour=:real)
	corr4 = QuAPI.Δt(bath0; N=4, t=0.4)
	@test_throws DimensionMismatch HybridizationKernel(lat2, corr4)

	# 虚时只有 τ 块
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernelτ = hybridization(latτ, bath0)
	@test size(kernelτ.blocks[(:τ, :τ)]) == (4, 4)

	# 混合九块
	latm = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	kernelm = hybridization(latm, bath0)
	for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
		@test haskey(kernelm.blocks, (b1, b2))
	end
end

@testset "API: hybriddynamics（IF 链）" begin
	lat = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real)
	for b in 1:2
		I = hybriddynamics(lat, bath0, PartialIF(); band=b, coupling_norm=1.0)
		@test I isa GM.GMPS && length(I) == length(lat) ÷ 2
		@test check_charge(I)
	end
	# band 越界显式拒绝
	@test_throws BoundsError hybriddynamics(lat, bath0, PartialIF(); band=3)

	# 虚时：PartialIF 与 ExactTTI 同一泛函的两种物化 → 同 Z
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	Iτ = hybriddynamics(latτ, bath0, PartialIF(); band=1)
	@test length(Iτ) == length(latτ) ÷ 2
	Ie = hybriddynamics(latτ, bath0, ExactTTI(); band=1)
	@test length(Ie) == length(latτ) ÷ 2
	@test integrate(Iτ) ≈ integrate(Ie) rtol = 1.0e-8
end

@testset "API: fit_exponentials 薄壳（ExpExp sweep 委托，plan kernel_and_fit 案 B）" begin
	# 精确指数和重建（衰减 + 复相位极点）；:left 约定 weights[i]·nodes[i]^(d−1)
	y = ComplexF64[0.7 * (0.5cis(0.3))^(d-1) + 0.3 * 0.2^(d-1) for d in 1:16]
	f = GTEMPO.fit_exponentials(y)
	@test length(f.nodes) == 2 && f.anchors == [:left, :left]
	@test all(abs(n) <= 1 for n in f.nodes) && f.relative_error <= 1e-10
	rec = [sum(f.weights[i] * f.nodes[i]^(d-1) for i in eachindex(f.nodes)) for d in 1:16]
	@test maximum(abs.(rec .- y)) / norm(y) <= 1e-10
	# 零序列 / L=1 早退（IF 域策略层）
	f0 = GTEMPO.fit_exponentials(zeros(ComplexF64, 5))
	@test isempty(f0.nodes) && isempty(f0.weights) && f0.relative_error == 0.0
	f1 = GTEMPO.fit_exponentials(ComplexF64[3.0])
	@test f1.nodes == [0.0] && f1.weights == [3.0] && f1.anchors == [:left]
	# 退化窗回退：线性序列铅笔不达 → 单位根精确插值（L ≤ maxterms）
	fr = GTEMPO.fit_exponentials(ComplexF64.(collect(1.0:8.0)))
	@test length(fr.nodes) == 8 && all(a === :left for a in fr.anchors)
	# 无回退资格（maxterms < L）→ sweep 的 ArgumentError 透传（带最优残差）
	@test_throws ArgumentError GTEMPO.fit_exponentials(ComplexF64.(collect(1.0:12.0)); maxterms=2)
	# rightanchors kwarg 已删（2026-11-17 恒双端锚定裁决）
	@test_throws MethodError GTEMPO.fit_exponentials(y; rightanchors=false)
	# 超临界极点折叠 :right（恒等式精确，节点 ≤ 1）
	yg = ComplexF64[1.2^(d-1) for d in 1:10]
	fg = GTEMPO.fit_exponentials(yg)
	@test fg.anchors == [:right] && abs(fg.nodes[1]) <= 1
	@test maximum(abs.([fg.weights[1] * fg.nodes[1]^(10-d) for d in 1:10] .- yg)) / norm(yg) <= 1e-10
end
