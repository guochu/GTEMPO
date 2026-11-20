# BCS 连续 bath 测试（test/bcs/continuousbath）：
#   ① 核结构：Δ=0 连续核 ≡ normal 连续核（连续积分噪声档 ~1e-11）；
#     tti 因子 fit 模式 smoke（Keldysh/热配对审计零残差——连续核同走 TTI 门）。
#   ② 实时演化 vs **IMB BCS-Toulouse 严格解**（toulouse_neq_greater_lesser，
#     Bogoliubov 机器）：TN 与锚用同一 discretebath(δw) 离散浴 ⟹ δw 误差
#     两侧抵消，容差只含 δt 离散化项（order=2 / order=1 分档）。
#   ③ Z 恒等式：U=0 线性耦合下 Keldysh 轨迹 Z=1 精确成立（order=2，δt² 档）。
println("---------- BCS continuous bath ----------")
using LinearAlgebra: I, dot

@testset "BCS: 连续核结构（Δ=0 ≡ normal + TTI 审计）" begin
	spec = spectrum_func()
	β = 1.0
	lat = GrassmannLattice(N=6, δt=0.05, bands=2, contour=:real, symmetry=Z2Symmetry())
	bath0 = bcsbath(fermionicbath(spec, β=β), Δ=0.0)
	kb0 = correlationfunction(bath0, lat)
	lat1 = GrassmannLattice(N=6, δt=0.05, bands=1, contour=:real)
	kn = correlationfunction(fermionicbath(spec, β=β), lat1)
	dev = 0.0
	for b1 in (:+, :-), b2 in (:+, :-)
		mk = kb0.blocks[(ContourBlockIndex(b1, 1, true), ContourBlockIndex(b2, 1, false))]
		mn = kn.blocks[(ContourBlockIndex(b1, 1, true), ContourBlockIndex(b2, 1, false))]
		dev = max(dev, maximum(abs.(Matrix(mk) - Matrix(mn))))
	end
	@test dev < 1e-10                     # 实测 5e-12（连续积分噪声档）
	@test kb0.anomalous == true           # 反常味键结构性存在（构造器全味键表）
	# 反常条目值恒零（Δ=0 ⟹ uv 谱权重 = 0）
	let deva = 0.0
		for (ka, ca, kb) in ((1, true, 2), (2, false, 1)), b1 in (:+, :-), b2 in (:+, :-)
			mk = kb0.blocks[(ContourBlockIndex(b1, ka, ca), ContourBlockIndex(b2, kb, ca))]
			deva = max(deva, maximum(abs.(Matrix(mk))))
		end
		@test deva == 0.0
	end
	# Δ=0.6：反常键存在、tti 因子全 fit 且审计零残差
	bath6 = bcsbath(fermionicbath(spec, β=β), Δ=0.6)
	kb6 = correlationfunction(bath6, lat)
	@test kb6.anomalous == true
	(; factors, audit) = tti_factors(kb6, lat)
	@test length(factors) == 16           # 4 味对 × 4 支块（实时）
	@test audit.keldysh_residual < 1e-12
	@test all(e.mode == :fit for e in audit.real_expansions)
	# 虚时：热配对审计 + 环族 fit
	lati = GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
	ki = correlationfunction(bath6, lati)
	fi = tti_factors(ki, lati)
	@test fi.audit.thermal_residual < 1e-12
	@test length(fi.factors) == 4         # 4 味对 × 环族 ESS 池
	@test all(e.mode == :fit for e in fi.audit.imag_expansions)
end

@testset "BCS: 连续 bath 实时演化 vs 同模型 BdG 锚（U=0）" begin
	β = 1.0; Δ = 0.6; ϵ_d = 0.7; T = 0.4; Nt = 8
	spec = spectrum_func()
	bathc = bcsbath(fermionicbath(spec, β=β), Δ=Δ)
	bathd = discretebath(bathc; δw=0.1)             # IMB 离散化（TN 与锚同一浴）
	model_ed = Toulouse(bathd; ϵ_d=ϵ_d)             # BCSToulouse（IMB Bogoliubov 严格解）
	model = AndersonIM(U=0.0, ϵ_d=ϵ_d)
	ρ0 = thermal_impurity_state(model, β)
	δt = T / Nt
	ts = (1:Nt) .* δt
	nsys = 1 / (exp(β * ϵ_d) + 1)                   # 杂质每带占据（U=0 闭式）
	anchor_gt, anchor_lt = toulouse_neq_greater_lesser(model_ed, ts; nsys=nsys)
	lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
	kernel = correlationfunction(bathd, lat)
	Is, _ = hybriddynamics(lat, kernel, PartialIF())
	Is = Is isa Tuple ? Is : (Is,)
	for order in (2, 1)
		K = order == 2 ?
			sysdynamics(lat, model; order=2, coupling_norm=2coupling_norm(bathc),
			            initial_state=ρ0) :
			sysdynamics(lat, model; initial_state=ρ0)
		Z = integrate(K, Is...)
		for i in 1:Nt
			g = greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=order,
			            coupling_norm=2coupling_norm(bathc), initial_state=ρ0)
			l = lesser(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=order,
			           coupling_norm=2coupling_norm(bathc), initial_state=ρ0)
			@test g ≈ anchor_gt[i] atol = (order == 2 ? 1e-3 : 5e-3)
			@test l ≈ anchor_lt[i] atol = (order == 2 ? 1e-3 : 5e-3)
		end
	end
end

@testset "BCS: Z 恒等式（U=0 线性耦合，连续 BCS bath）" begin
	β = 1.0; Δ = 0.6; ϵ_d = 0.7
	spec = spectrum_func()
	bathc = bcsbath(fermionicbath(spec, β=β), Δ=Δ)
	bathd = discretebath(bathc; δw=0.1)             # IMB 离散化（TN 与锚同一浴）
	model = AndersonIM(U=0.0, ϵ_d=ϵ_d)
	ρ0 = thermal_impurity_state(model, β)
	for (Nt, atolz) in ((4, 2e-2), (8, 5e-3))
		δt = 0.4 / Nt
		lat = GrassmannLattice(N=Nt, δt=δt, bands=2, contour=:real, symmetry=Z2Symmetry())
		kernel = correlationfunction(bathd, lat)
		Is, _ = hybriddynamics(lat, kernel, PartialIF())
		Is = Is isa Tuple ? Is : (Is,)
		K = sysdynamics(lat, model; order=2, coupling_norm=2coupling_norm(bathc),
		                initial_state=ρ0)
		# 实测 5.4e-3 / 1.4e-3（连续谱 δt² 档，常数大于离散谱）
		@test abs(integrate(K, Is...) - 1) < atolz
	end
end
