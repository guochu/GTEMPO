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
		blk = kernel.blocks[(ContourBlockIndex(b1, 1, true), ContourBlockIndex(b2, 1, false))]
		@test size(blk) == (4, 4)
		@test all(isfinite, blk)
	end
	# BCS/normal 同类型实证（F5）：同一 struct、不同构造方法
	bath_bcs = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=β0), Δ=0.3)
	lat_b = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real, symmetry=Z2Symmetry())
	kbcs = correlationfunction(bath_bcs, lat_b)
	@test kbcs isa HybridizationKernel && kbcs.bands == [1, 2]
	@test kbcs.coupling_norm == 2 * coupling_norm(bath_bcs)   # 带对各 α → 2×
	@test haskey(kbcs.blocks, (ContourBlockIndex(:+, 1, true), ContourBlockIndex(:-, 2, true)))  # ud 反常键（双腿 ā）
	@test haskey(kbcs.blocks, (ContourBlockIndex(:-, 2, false), ContourBlockIndex(:+, 1, false)))  # du 反常键（双腿 a）
	# orbital=nothing：全部 Nambu 对共享同一套块对象（零拷贝别名，normal 多带同款）
	bath_bcs4 = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=β0), Δ=0.3)
	lat_b4 = GrassmannLattice(N=4, δt=0.1, bands=4, contour=:real, symmetry=Z2Symmetry())
	kb_all = correlationfunction(bath_bcs4, lat_b4)
	@test kb_all.bands == [1, 2, 3, 4]
	@test kb_all.coupling_norm == 4 * coupling_norm(bath_bcs4)
	for (dba, dbb) in ((1, 1), (1, 2), (2, 2), (2, 1))
		ca, cb = GTEMPO._bcs_leg_conjs(dba, dbb)
		for b1 in (:+, :-), b2 in (:+, :-)
			@test kb_all.blocks[(ContourBlockIndex(b1, 2 + dba, ca), ContourBlockIndex(b2, 2 + dbb, cb))] ===
			      kb_all.blocks[(ContourBlockIndex(b1, dba, ca), ContourBlockIndex(b2, dbb, cb))]
		end
	end
	kb_o2 = correlationfunction(bath_bcs4, lat_b4; orbital=2)
	@test kb_o2.bands == [3, 4] && kb_o2.coupling_norm == 2 * coupling_norm(bath_bcs4)
	@test Matrix(kb_o2.blocks[(ContourBlockIndex(:-, 3, true), ContourBlockIndex(:+, 4, true))]) ≈
	      Matrix(kb_all.blocks[(ContourBlockIndex(:-, 1, true), ContourBlockIndex(:+, 2, true))]) atol = 0.0

	# 虚时只有 τ 块
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernelτ = correlationfunction(bath0, latτ)
	@test size(kernelτ.blocks[(ContourBlockIndex(:τ, 1, true), ContourBlockIndex(:τ, 1, false))]) == (4, 4)

	# 混合九块 + τ↔rt 构造期折叠（存储 = 指数系数终值）
	latm = GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.25, contour=:mixed)
	kernelm = correlationfunction(bath0, latm)
	for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
		@test haskey(kernelm.blocks, (ContourBlockIndex(b1, 1, true), ContourBlockIndex(b2, 1, false)))
	end
end

@testset "API: changebands（带重排）" begin
	lat = GrassmannLattice(N=4, δt=0.1, contour=:real)
	k1 = correlationfunction(bath0, lat; band=1)
	k12 = changebands(k1, [1 => 2])
	@test k12.bands == [2]
	@test k12.blocks[(ContourBlockIndex(:+, 2, true), ContourBlockIndex(:+, 2, false))] ===
	      k1.blocks[(ContourBlockIndex(:+, 1, true), ContourBlockIndex(:+, 1, false))]   # 零改动（对象别名）
	@test k12.coupling_norm == k1.coupling_norm
	# 单带核扩多带（免重建 corr 的多带装配路线）
	k1v2 = changebands(k1, [1 => 1, 1 => 2])
	@test k1v2.bands == [1, 2]
	@test k1v2.blocks[(ContourBlockIndex(:-, 2, true), ContourBlockIndex(:-, 2, false))] ===
	      k1v2.blocks[(ContourBlockIndex(:-, 1, true), ContourBlockIndex(:-, 1, false))]
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
	Ma = Matrix(ka.blocks[(ContourBlockIndex(:+, 1, true), ContourBlockIndex(:-, 1, false))])
	Mb = Matrix(kb.blocks[(ContourBlockIndex(:+, 1, true), ContourBlockIndex(:-, 1, false))])
	@test Matrix(kab.blocks[(ContourBlockIndex(:+, 1, true), ContourBlockIndex(:-, 1, false))]) ≈ Ma + Mb atol = 0.0
	# 缺键透传（单带 vs 双带核相加）
	kb_multi = changebands(correlationfunction(b2, lat), [1 => 2])
	ksum = ka + kb_multi
	@test ksum.bands == [1, 2]
	@test ksum.blocks[(ContourBlockIndex(:+, 2, true), ContourBlockIndex(:-, 2, false))] ===
	      kb_multi.blocks[(ContourBlockIndex(:+, 2, true), ContourBlockIndex(:-, 2, false))]
end

@testset "API: hybriddynamics（IF 链）" begin
	lat = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real)
	kernel = correlationfunction(bath0, lat)
	Is, _ = hybriddynamics(lat, kernel, PartialIF())
	@test Is isa Tuple && length(Is) == 2                 # 多带 kernel → 逐带 IF 序
	for I in Is
		@test I isa GM.GMPS && length(I) == length(lat) ÷ 2
		@test check_charge(I)
	end
	# 原地版：多带核逐带从 copy(gmps) 起步（多带返回元组，gmps 本身不变）
	Ik, diagk = hybriddynamics!(vacuumstate(lat), lat, kernel, PartialIF())
	@test Ik isa Tuple && length(Ik) == length(Is)
	@test all(i -> isapprox(integrate(Ik[i]), integrate(Is[i]); rtol=1e-12), eachindex(Is))
	@test diagk[:truncerrs] isa Vector{Float64}
	@test haskey(diagk, :algorithm) && haskey(diagk, :truncation)
	# 单带 kernel → 单 IF + 原地恒等 + diagnostics 第二输出
	k2 = correlationfunction(bath0, lat; band=2)
	g2 = vacuumstate(lat)
	I2, diag2 = hybriddynamics!(g2, lat, k2, PartialIF())
	@test I2 === g2
	@test diag2[:truncerrs] isa Vector{Float64}
	# band 越界显式拒绝
	@test_throws BoundsError correlationfunction(bath0, lat; band=3)

	# 虚时：PartialIF 与 ExactTTI 同一泛函的两种物化 → 同 Z
	latτ = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
	kernelτ = correlationfunction(bath0, latτ)
	Iτ, _ = hybriddynamics(latτ, kernelτ, PartialIF())
	@test length(Iτ) == length(latτ) ÷ 2
	@test TK.scalartype(Iτ) == Float64     # 虚时实谱链 = Float64（lattice Float64 × kernel Float64）
	Ie, _ = hybriddynamics(latτ, kernelτ, ExactTTI())
	@test length(Ie) == length(latτ) ÷ 2
	@test TK.scalartype(Ie) == ComplexF64  # ExactTTI 拟合层恒复（实链装配入口复化）
	@test integrate(Iτ) ≈ integrate(Ie) rtol = 1.0e-8
end

@testset "API: fit_exponentials（alg_expan 面：AnchoredMatrixPencil + 有界回退）" begin
	# 精确指数和重建（衰减 + 复相位极点）；y[d] = Σ w·q^(d−1)
	D = GTEMPO.Defaults
	y = ComplexF64[0.7 * (0.5cis(0.3))^(d-1) + 0.3 * 0.2^(d-1) for d in 1:16]
	f = GTEMPO.fit_exponentials(D.alg_expan(), y)
	@test length(f.nodes) == 2 && all(abs(n) <= 1 for n in f.nodes) && f.relative_error <= 1e-10
	rec = [sum(f.weights[i] * f.nodes[i]^(d-1) for i in eachindex(f.nodes)) for d in 1:16]
	@test maximum(abs.(rec .- y)) / norm(y) <= 1e-10
	# 零序列 / L=1 早退（IF 域策略层）
	f0 = GTEMPO.fit_exponentials(D.alg_expan(), zeros(ComplexF64, 5))
	@test isempty(f0.nodes) && isempty(f0.weights) && f0.relative_error == 0.0
	f1 = GTEMPO.fit_exponentials(D.alg_expan(), ComplexF64[3.0])
	@test f1.nodes == [0.0] && f1.weights == [3.0] && f1.relative_error == 0.0
	# 契约断言：err ≤ alg.tol（铅笔收敛或单位根回退，节点数随收敛路径 2–8）
	fr = GTEMPO.fit_exponentials(D.alg_expan(), ComplexF64.(collect(1.0:8.0)))
	@test fr.relative_error <= 1e-10 && 2 <= length(fr.nodes) <= 8
	# 线性 50 点：铅笔残差 1.03e-10 边缘超 tol 且 L > alg.n 无回退 → 抛
	@test_throws ArgumentError GTEMPO.fit_exponentials(D.alg_expan(), ComplexF64.(collect(1.0:50.0)))
	# 无规律噪声序列：8 点复序列恰被 4 极点精确表示（16 dof = 16 约束）——
	# 铅笔机器精度收敛（超临界节点直出、代数精确），回退不触发（无需求即正确）
	using Random
	rng = Random.MersenneTwister(42)
	yv = ComplexF64.(randn(rng, 8) .+ im .* randn(rng, 8))
	fv = GTEMPO.fit_exponentials(D.alg_expan(), yv)
	@test fv.relative_error <= 1e-12
	# 超临界极点：节点直出（|q| = 1.2 未折叠），y[d] = Σ w·q^(d−1) 重建精确——
	# 下游 Pairing/ESS 装配对任意模长节点代数精确（守卫已随锚点概念退役）
	yg = ComplexF64[1.2^(d-1) for d in 1:10]
	fg = GTEMPO.fit_exponentials(D.alg_expan(), yg)
	@test length(fg.nodes) == 1 && abs(fg.nodes[1] - 1.2) <= 1e-10
	@test maximum(abs.([fg.weights[1] * fg.nodes[1]^(d-1) for d in 1:10] .- yg)) / norm(yg) <= 1e-10
end
