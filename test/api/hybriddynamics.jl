using LinearAlgebra: norm
import TensorKit as TK
β0 = 1.0
bath0 = discretefermionicbath([1.0], [0.5]; β=β0)

@testset "API: HybridizationKernel（W2 重定义）" begin
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	kernel = correlationfunction(bath0, lat)
	@test kernel isa HybridizationKernel
	@test TK.scalartype(typeof(kernel)) == ComplexF64
	@test TK.scalartype(kernel) == ComplexF64          # 泛型回退覆盖实例（仅 Type 定义）
	@test kernel.bands == [1]
	@test kernel.coupling_norm == coupling_norm(bath0)  # 单带 = bath α
	for b1 in (:+, :-), b2 in (:+, :-)
		blk = kernel.blocks[(ContourBlockIndex(b1, 1), ContourBlockIndex(b2, 1))]
		@test size(blk) == (4, 4)
		@test all(isfinite, blk)
	end
	# BCS/normal 同类型实证（F5）：同一 struct、不同构造方法
	bath_bcs = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=β0), Δ=0.3)
	lat_b = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real, symmetry=Z2Symmetry())
	kbcs = correlationfunction(bath_bcs, lat_b)
	@test kbcs isa HybridizationKernel && kbcs.bands == [1, 2]
	@test kbcs.coupling_norm == 2 * coupling_norm(bath_bcs)   # 带对各 α → 2×
	@test haskey(kbcs.blocks, (ContourBlockIndex(:+, 1), ContourBlockIndex(:-, 2)))  # ud 反常键
	@test haskey(kbcs.blocks, (ContourBlockIndex(:-, 2), ContourBlockIndex(:+, 1)))  # du 反常键

	# 虚时只有 τ 块
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernelτ = correlationfunction(bath0, latτ)
	@test size(kernelτ.blocks[(ContourBlockIndex(:τ, 1), ContourBlockIndex(:τ, 1))]) == (4, 4)

	# 混合九块 + τ↔rt 构造期折叠（存储 = 指数系数终值）
	latm = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	kernelm = correlationfunction(bath0, latm)
	for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
		@test haskey(kernelm.blocks, (ContourBlockIndex(b1, 1), ContourBlockIndex(b2, 1)))
	end
end

@testset "API: changebands（带重排）" begin
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	k1 = correlationfunction(bath0, lat; band=1)
	k12 = changebands(k1, [1 => 2])
	@test k12.bands == [2]
	@test k12.blocks[(ContourBlockIndex(:+, 2), ContourBlockIndex(:+, 2))] ===
	      k1.blocks[(ContourBlockIndex(:+, 1), ContourBlockIndex(:+, 1))]   # 零改动（对象别名）
	@test k12.coupling_norm == k1.coupling_norm
	# 单带核扩多带（免重建 corr 的多带装配路线）
	k1v2 = changebands(k1, [1 => 1, 1 => 2])
	@test k1v2.bands == [1, 2]
	@test k1v2.blocks[(ContourBlockIndex(:-, 2), ContourBlockIndex(:-, 2))] ===
	      k1v2.blocks[(ContourBlockIndex(:-, 1), ContourBlockIndex(:-, 1))]
	@test_throws ArgumentError changebands(k1, [3 => 2])          # 源带不在核
	@test_throws ArgumentError changebands(k1, [1 => 1, 1 => 1])  # 目标带冲突
end

@testset "API: kernel 相加（多 bath 耦合）" begin
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	b1 = discretefermionicbath([1.0], [0.5]; β=1.0)
	b2 = discretefermionicbath([2.0], [0.3]; β=1.0)
	ka = correlationfunction(b1, lat)
	kb = correlationfunction(b2, lat)
	kab = ka + kb
	@test kab.coupling_norm ≈ coupling_norm(b1) + coupling_norm(b2)
	@test kab.bands == [1]
	# CorrelationMatrix 值的 Toeplitz 向量相加（逐位 = 稠密和）
	Ma = Matrix(ka.blocks[(ContourBlockIndex(:+, 1), ContourBlockIndex(:-, 1))])
	Mb = Matrix(kb.blocks[(ContourBlockIndex(:+, 1), ContourBlockIndex(:-, 1))])
	@test Matrix(kab.blocks[(ContourBlockIndex(:+, 1), ContourBlockIndex(:-, 1))]) ≈ Ma + Mb atol = 0.0
	# 缺键透传（单带 vs 双带核相加）
	kb_multi = changebands(correlationfunction(b2, lat), [1 => 2])
	ksum = ka + kb_multi
	@test ksum.bands == [1, 2]
	@test ksum.blocks[(ContourBlockIndex(:+, 2), ContourBlockIndex(:-, 2))] ===
	      kb_multi.blocks[(ContourBlockIndex(:+, 2), ContourBlockIndex(:-, 2))]
end

@testset "API: hybriddynamics（IF 链）" begin
	lat = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real)
	kernel = correlationfunction(bath0, lat)
	Is = hybriddynamics(lat, kernel, PartialIF())
	@test Is isa Tuple && length(Is) == 2                 # 多带 kernel → 逐带 IF 序
	for I in Is
		@test I isa GM.GMPS && length(I) == length(lat) ÷ 2
		@test check_charge(I)
	end
	Ik = hybriddynamics(lat, kernel; alg=PartialIF())          # 关键字形委托
	@test Ik isa Tuple && length(Ik) == length(Is)
	@test all(i -> isapprox(integrate(Ik[i]), integrate(Is[i]); rtol=1e-12), eachindex(Is))
	# 单带 kernel → 单 IF
	k2 = correlationfunction(bath0, lat; band=2)
	I2 = hybriddynamics(lat, k2, PartialIF())
	@test I2 isa GM.GMPS
	# band 越界显式拒绝
	@test_throws BoundsError correlationfunction(bath0, lat; band=3)

	# 虚时：PartialIF 与 ExactTTI 同一泛函的两种物化 → 同 Z
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernelτ = correlationfunction(bath0, latτ)
	Iτ = hybriddynamics(latτ, kernelτ, PartialIF())
	@test length(Iτ) == length(latτ) ÷ 2
	Ie = hybriddynamics(latτ, kernelτ, ExactTTI())
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
