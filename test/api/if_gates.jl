# IF 装配闸门测试（plan/if_expmult_backend_plan.md Phase 2 验收，2026-11-18
# 旧装配档案退役后为互换/守卫闸门；2026-11-19 PolarIF 对拍并入各算例——
# 同一 PartialIF 装配态复用，免重复装配）：
#   ① 互换对拍：ExactTTI ≡ PartialIF（虚时/混合机器精度——混合格 Pairing 直读
#     路线无拟合；实时 = 拟合残差档）+ PolarIF ≡ PartialIF（同一 fullfactor
#     第二种精确因子化，差 = 压缩路径噪声档，逐例实测锚）；
#   ② 跨对称性互换（U1 normal：fullfactor 仅 ηη̄ 条目可物化）；
#   ③ 守卫：实时增长极点响抛、BCS × ExactTTI 拒绝（D7）、band 范围守卫；
#   ED 对拍由既有套件承接（normalbath/bcs/multiorbital + exactttiif golden）。
@testset "IF 装配闸门（W3/W4 后端统一）" begin
	gd(x, y) = maximum(abs.(GM.Dense.todense(x).data - GM.Dense.todense(y).data))

	# ---- ① 虚时互换（单模浴，μ=0 / μ≠0，拟合精确 → 机器精度）----
	for (bath, label) in ((discretefermionicbath([1.0], [0.5]; β=1.0), "μ=0"),
	                      (fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0, μ=0.7), "μ=0.7"))
		lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)   # N=3 触发热配对审计轮盘（已知）
		kernel = correlationfunction(bath, lat)
		full, _ = hybriddynamics(lat, kernel, PartialIF())
		exact, _ = hybriddynamics(lat, kernel, ExactTTI())
		polar, _ = hybriddynamics(lat, kernel, PolarIF())
		@test gd(exact, full) < 1e-12
		@test gd(polar, full) < 1e-13     # 实测 ~1e-15
	end

	# ---- ① 实时互换（连续谱阻尼核；缺省 alg_expan = OverDeterminedProny(15,1e-4)
	#       对短核列残差 ~1e-6 档 ⟹ 互换 ~2.5e-3；高精度配
	#       ExactTTI(alg_expan = Defaults.alg_expan(tol=1e-10, n=40)) 或锚定铅笔）----
	let lat = GrassmannLattice(N=3, δt=0.15, contour=:real)
		spec_c = spectrum(ϵ -> sqrt(1 - (ϵ / 10)^2) / π, lb=-10, ub=10)
		bath = fermionicbath(spec_c, β=1.0)
		kernel = correlationfunction(bath, lat)
		exact, _ = hybriddynamics(lat, kernel, ExactTTI())
		rowchain, _ = hybriddynamics(lat, kernel, PartialIF())
		polar, _ = hybriddynamics(lat, kernel, PolarIF())   # raw 直读无拟合档
		@test gd(exact, rowchain) < 1e-2
		@test gd(polar, rowchain) < 1e-13                   # 实测 1.2e-15（缩格后压缩无损）
	end

	# ---- ③ DiracDelta 实时边界极点（|q|=1）可装配（2026-11-18 alg_expan 语义：
	#       无锚点折叠 ⟹ 旧 :right 守卫退役；单模核 Prony 精确 ⟹ 机器精度互换）----
	let lat = GrassmannLattice(N=3, δt=0.15, contour=:real)
		bath = fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0)
		kernel = correlationfunction(bath, lat)
		exact, _ = hybriddynamics(lat, kernel, ExactTTI())
		rowchain, _ = hybriddynamics(lat, kernel, PartialIF())
		polar, _ = hybriddynamics(lat, kernel, PolarIF())
		@test gd(exact, rowchain) < 1e-12
		@test gd(polar, rowchain) < 1e-13
	end

	# ---- ① 混合互换（Pairing 直读路线——ExactTTI 机器精度；PolarIF 压缩噪声档）----
	let latm = GrassmannLattice(Nt=2, δt=0.15, Nτ=3, δτ=0.25, contour=:mixed)
		spec_c = spectrum(ϵ -> sqrt(1 - (ϵ / 10)^2) / π, lb=-10, ub=10)
		bath = fermionicbath(spec_c, β=latm.β)
		kernel = correlationfunction(bath, latm)
		exact, _ = hybriddynamics(latm, kernel, ExactTTI())
		rowchain, _ = hybriddynamics(latm, kernel, PartialIF())
		polar, _ = hybriddynamics(latm, kernel, PolarIF())
		@test gd(exact, rowchain) < 1e-12
		@test gd(polar, rowchain) < 1e-13   # 实测 5.3e-15
	end

	# ---- ② 跨对称性互换（U1 normal：fullfactor 仅 ηη̄ 条目可物化）----
	let lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag, symmetry=U1Symmetry())
		bath = discretefermionicbath([1.0], [0.5]; β=1.0)
		kernel = correlationfunction(bath, lat)
		exact, _ = hybriddynamics(lat, kernel, ExactTTI())
		rowchain, _ = hybriddynamics(lat, kernel, PartialIF())
		polar, _ = hybriddynamics(lat, kernel, PolarIF())
		@test gd(exact, rowchain) < 1e-12
		@test gd(polar, rowchain) < 1e-13
	end

	# ---- ③ BCS × ExactTTI 拒绝（D7）+ 多带核 fullfactor/tti_* 范围守卫 ----
	let latb = GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
		bathb = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0), Δ=0.4)
		kernelb = correlationfunction(bathb, latb)
	@test_throws ArgumentError hybriddynamics(latb, kernelb, ExactTTI())
	@test fullfactor(kernelb, latb) isa GM.AnomalousPairing      # BCS 轨道核（anomalous）OK
	bathn = discretefermionicbath([1.0], [0.5]; β=1.0)
	kfull = correlationfunction(bathn, latb)                     # 多带 normal 核
	ff = fullfactor(kfull, latb)                                 # 泛化：normal = Pairing 单因子全带
	@test ff isa GM.Pairing && ff.asites == ff.bsites && length(ff.asites) == length(latb) ÷ 2
	@test_throws ArgumentError tti_imag_factors(kfull, latb)     # tti_* 仍限单带核
	end
end

@testset "IF 装配闸门：PolarIF（BCS 反常 + dense oracle + rtol smoke）" begin
	# plan/polarif_plan.md T1/T2（对拍算例已并入 W3/W4 闸门复用装配态）：
	# 此处只留 W3/W4 未覆盖的面——AnomalousPairing 味扇区折叠（BCS）、跨带随机核
	# dense oracle（vs 精确 exp）、恒单链契约、rtol 截秩 smoke。
	gd(x, y) = maximum(abs.(GM.Dense.todense(x).data - GM.Dense.todense(y).data))

	# ---- BCS 虚时反常（AnomalousPairing 味扇区折叠路径；实测 ~1e-8）----
	let latb = GrassmannLattice(N=3, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
		bathb = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=latb.β), Δ=0.4)
		kernelb = correlationfunction(bathb, latb)
		full, _ = hybriddynamics(latb, kernelb, PartialIF())
		polar, _ = hybriddynamics(latb, kernelb, PolarIF())
		@test gd(polar, full) < 1e-13       # 实测 1.2e-14（N=3 压缩无损；β 同源免告警）
	end

	# ---- 跨带随机核 dense oracle（vs 精确 exp；实测 ≤ 3.2e-12）+ 恒单链契约 ----
	DT = GM.Dense
	import Random
	function rand_kernel(rng, lat, bands_list)
		T = ComplexF64
		blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}()
		brs = branches(lat)
		for b1 in bands_list, b2 in bands_list, br1 in brs, br2 in brs
			n1 = br1 === :τ ? lat.kτ : lat.kt
			n2 = br2 === :τ ? lat.kτ : lat.kt
			blocks[(ContourBlockIndex(br1, b1, true), ContourBlockIndex(br2, b2, false))] =
				randn(rng, T, n1, n2)
		end
		HybridizationKernel(blocks, 0.0, sort(unique(bands_list)))
	end
	gd2(x, y) = maximum(abs.(x .- y))
	# ⚠️ dense oracle 按 2^(2·站数) 分配——案例站数上限 8；D=512 覆盖 M ≤ 8 无截断。
	alg = PolarIF(trunc=GTEMPO.Defaults.alg_trunc_if(D=512))
	for (B, kbands, N) in ((2, [1, 2], 3), (3, [1, 2, 3], 2))
		rng = Random.MersenneTwister(hash((B, tuple(kbands...), N)))
		lat = GrassmannLattice(N=N, δτ=0.25, bands=B, contour=:imag)
		kernel = rand_kernel(rng, lat, kbands)
		state, diag = hybriddynamics!(vacuumstate(lat), lat, kernel, alg)
		oracle = DT.tositewise(DT.exp(DT.todense(fullfactor(kernel, lat), length(lat) ÷ 2)))
		@test gd2(DT.tositewise(GM.Dense.todense(state)), oracle) < 1e-10
	end
	let rng = Random.MersenneTwister(777)
		lat = GrassmannLattice(N=2, δt=0.15, bands=2, contour=:real)
		kernel = rand_kernel(rng, lat, [1, 2])
		state, _ = hybriddynamics!(vacuumstate(lat), lat, kernel, alg)
		oracle = DT.tositewise(DT.exp(DT.todense(fullfactor(kernel, lat), length(lat) ÷ 2)))
		@test gd2(DT.tositewise(GM.Dense.todense(state)), oracle) < 1e-10
	end
	# 恒单链契约（继承 PartialIF）：多带核 PolarIF 原地恒等返回单 GMPS
	let rng = Random.MersenneTwister(555)
		lat = GrassmannLattice(N=3, δτ=0.25, bands=2, contour=:imag)
		kernel = rand_kernel(rng, lat, [1, 2])
		g = vacuumstate(lat)
		state, _ = hybriddynamics!(g, lat, kernel, PolarIF())
		@test state === g && state isa GM.GMPS && length(state) == length(lat) ÷ 2
	end

	# ---- T2 rtol 截秩 smoke（数值锚实测回填：rtol=0.1 不掉秩=逐位同 rtol=0；
	#       rtol=0.5 掉半秩质量 ⟹ 可观测差 4.0e-3）----
	let lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
		bath = discretefermionicbath([1.0], [0.5]; β=1.0)
		kernel = correlationfunction(bath, lat)
		exact, _ = hybriddynamics(lat, kernel, PolarIF())
		quiet, _ = hybriddynamics(lat, kernel, PolarIF(rtol=0.1))
		@test gd(quiet, exact) == 0.0                     # 秩谱尾 < 0.1 → 无截 → 精确同
		truncd, diag = hybriddynamics(lat, kernel, PolarIF(rtol=0.5))
		d = gd(truncd, exact)
		@test 0 < d < 1e-1                                # 实测 4.0e-3
		@test isfinite(diag[:truncerr])
		@test isfinite(integrate(truncd))
	end
end

@testset "泛化 kernel：随机跨带耦合 PartialIF（dense oracle）" begin
	# 一般 HybridizationKernel 可含跨带耦合（Δ_b1b2 ≠ 0）——二次型不可按带分解，
	# hybriddynamics 走单链整体装配（fullfactor 紧凑站点表）；正确性 = 装配态
	# vs 严格 exp(todense(fullfactor)) 的 GM dense oracle（DT.tositewise 契约）。
	DT = GM.Dense
	import Random
	function rand_kernel(rng, lat, bands_list)
		T = ComplexF64
		blocks = Dict{Tuple{ContourBlockIndex,ContourBlockIndex},Union{Matrix{T},CorrelationMatrix{T}}}()
		brs = branches(lat)
		for b1 in bands_list, b2 in bands_list, br1 in brs, br2 in brs
			n1 = br1 === :τ ? lat.kτ : lat.kt
			n2 = br2 === :τ ? lat.kτ : lat.kt
			blocks[(ContourBlockIndex(br1, b1, true), ContourBlockIndex(br2, b2, false))] =
				randn(rng, T, n1, n2)          # 一般稠密核（非 Toeplitz，Matrix 值）
		end
		HybridizationKernel(blocks, 0.0, sort(unique(bands_list)))
	end
	gd(x, y) = maximum(abs.(x .- y))

	# ⚠️ dense oracle 按 2^(2·站数) 分配——案例站数上限 8（2^16 = 65k 项），
	# 更大格走 ED/互换闸门验算，勿在此扩 N（内存爆炸点）。
	# 随机稠密二次型的精确指数中点键维上限 = 2^M——D=512 覆盖 M ≤ 8 案例的无截断
	# 验证（缺省 D=128 在 M ≥ 8 会截断，属算法语义非 bug）。
	alg = PartialIF(trunc=GTEMPO.Defaults.alg_trunc_if(D=512))
	for (B, kbands, N) in ((2, [1, 2], 3), (3, [1, 2, 3], 2), (4, [2, 3], 2), (4, [1, 3], 2))
		rng = Random.MersenneTwister(hash((B, tuple(kbands...), N)))
		lat = GrassmannLattice(N=N, δτ=0.25, bands=B, contour=:imag)
		kernel = rand_kernel(rng, lat, kbands)
		state, diag = hybriddynamics!(vacuumstate(lat), lat, kernel, alg)
		oracle = DT.tositewise(DT.exp(DT.todense(fullfactor(kernel, lat), length(lat) ÷ 2)))
		@test gd(DT.tositewise(GM.Dense.todense(state)), oracle) < 1e-10
	end

	# 实时轮廓（± 分支槽位数学）+ 跨 2 带全耦合
	let rng = Random.MersenneTwister(777)
		lat = GrassmannLattice(N=2, δt=0.15, bands=2, contour=:real)
		kernel = rand_kernel(rng, lat, [1, 2])
		state, _ = hybriddynamics!(vacuumstate(lat), lat, kernel, alg)
		oracle = DT.tositewise(DT.exp(DT.todense(fullfactor(kernel, lat), length(lat) ÷ 2)))
		@test gd(DT.tositewise(GM.Dense.todense(state)), oracle) < 1e-10
	end

	# 恒单链契约：带对角多带核与跨带核同为单个多带因子装配（逐带元组退役）
	let rng = Random.MersenneTwister(99)
		lat = GrassmannLattice(N=3, δτ=0.25, bands=3, contour=:imag)
		diag_kernel = rand_kernel(rng, lat, [1, 2, 3])
		filter!(p -> p[1][1].band == p[1][2].band, diag_kernel.blocks)   # 只留带对角键
		@test diag_kernel.anomalous == false && fullfactor(diag_kernel, lat) isa GM.Pairing
		Is, _ = hybriddynamics(lat, diag_kernel, PartialIF())
		@test Is isa GM.GMPS
		kernel = rand_kernel(rng, lat, [1, 2, 3])
		state, _ = hybriddynamics(lat, kernel, PartialIF())
		@test state isa GM.GMPS
	end
end
