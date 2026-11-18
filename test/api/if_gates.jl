# IF 装配闸门测试（plan/if_expmult_backend_plan.md Phase 2 验收，2026-11-18
# 旧装配档案退役后为互换/守卫闸门）：
#   ① 互换对拍：ExactTTI ≡ PartialIF（虚时/混合机器精度——混合格 Pairing 直读
#     路线无拟合；实时 = 拟合残差档，与旧路径同源同 rtol）；
#   ② 跨对称性互换（U1 normal：fullfactor 仅 ηη̄ 条目可物化）；
#   ③ 守卫：实时增长极点响抛、BCS × ExactTTI 拒绝（D7）、band 范围守卫；
#   ED 对拍由既有套件承接（normalbath/bcs/multiorbital + exactttiif golden）。
@testset "IF 装配闸门（W3/W4 后端统一）" begin
	gd(x, y) = maximum(abs.(GM.todense(x; order=:msb) - GM.todense(y; order=:msb)))

	# ---- ① 虚时互换（单模浴，μ=0 / μ≠0，拟合精确 → 机器精度）----
	for (bath, label) in ((discretefermionicbath([1.0], [0.5]; β=1.0), "μ=0"),
	                      (fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0, μ=0.7), "μ=0.7"))
		lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
		kernel = correlationfunction(bath, lat)
		full = influence_functional(lat, kernel; alg=PartialIF()).state
		exact = influence_functional(lat, kernel; alg=ExactTTI()).state
		@test gd(exact, full) < 1e-12
	end

	# ---- ① 实时互换（连续谱阻尼核；偏差 = 拟合残差档，旧路径同源同 rtol）----
	let lat = GrassmannLattice(N=4, δt=0.15, contour=:real)
		spec_c = spectrum(ϵ -> sqrt(1 - (ϵ / 10)^2) / π, lb=-10, ub=10)
		bath = fermionicbath(spec_c, β=1.0)
		kernel = correlationfunction(bath, lat)
		exact = influence_functional(lat, kernel; alg=ExactTTI()).state
		rowchain = influence_functional(lat, kernel; alg=PartialIF()).state
		@test gd(exact, rowchain) < 1e-9          # cross_rtol=1e-10 拟合残差档
	end

	# ---- ③ 实时增长极点守卫（DiracDelta |q|=1 边界极点响抛；折叠朝向随链长/
	#       位型抖动轮盘——N=3 为本机实测 :right 位型）----
	let lat = GrassmannLattice(N=3, δt=0.15, contour=:real)
		bath = fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0)
		kernel = correlationfunction(bath, lat)
		@test_throws ArgumentError influence_functional(lat, kernel; alg=ExactTTI())
	end

	# ---- ① 混合互换（Pairing 直读路线——机器精度）----
	let latm = GrassmannLattice(Nt=3, δt=0.15, Nτ=4, δτ=0.25, contour=:mixed)
		spec_c = spectrum(ϵ -> sqrt(1 - (ϵ / 10)^2) / π, lb=-10, ub=10)
		bath = fermionicbath(spec_c, β=latm.β)
		kernel = correlationfunction(bath, latm)
		exact = influence_functional(latm, kernel; alg=ExactTTI()).state
		rowchain = influence_functional(latm, kernel; alg=PartialIF()).state
		@test gd(exact, rowchain) < 1e-12
	end

	# ---- ② 跨对称性互换（U1 normal：fullfactor 仅 ηη̄ 条目可物化）----
	let lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag, symmetry=U1Symmetry())
		bath = discretefermionicbath([1.0], [0.5]; β=1.0)
		kernel = correlationfunction(bath, lat)
		exact = influence_functional(lat, kernel; alg=ExactTTI()).state
		rowchain = influence_functional(lat, kernel; alg=PartialIF()).state
		@test gd(exact, rowchain) < 1e-12
	end

	# ---- ③ BCS × ExactTTI 拒绝（D7）+ 多带核 fullfactor/tti_* 范围守卫 ----
	let latb = GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
		bathb = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0), Δ=0.4)
		kernelb = correlationfunction(bathb, latb)
		@test_throws ArgumentError influence_functional(latb, kernelb; alg=ExactTTI())
		@test fullfactor(kernelb, latb) isa GM.AnomalousPairing      # BCS 轨道核 OK
		bathn = discretefermionicbath([1.0], [0.5]; β=1.0)
		kfull = correlationfunction(bathn, latb)                     # 多带 normal 核
		@test_throws ArgumentError fullfactor(kfull, latb)
		@test_throws ArgumentError tti_imag_factors(kfull, latb)
	end
end
