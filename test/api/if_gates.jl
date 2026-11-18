# IF 装配闸门测试（plan/if_expmult_backend_plan.md Phase 2 验收）：
#   ① 每因子 dense 对拍：fullfactor/tti_* 产物 vs 现行（冻结档案）因子链
#     ——机器精度（拟合精确核）；
#   ② 互换对拍：ExactTTI ≡ PartialIF（虚时机器精度；实时 = 拟合残差档；
#     混合 = 机器精度——混合格 Pairing 路线无拟合）；
#   ③ ED 对拍由既有套件承接（normalbath/bcs/multiorbital + exactttiif fixture）。
@testset "IF 装配闸门（W3/W4 后端统一）" begin
	spec = TK.sectortype(Z2Symmetry())
	gd(x, y) = maximum(abs.(GM.todense(x; order=:msb) - GM.todense(y; order=:msb)))

	# ---- ① 虚时每因子 dense 对拍（单模浴，μ=0 / μ≠0，拟合精确）----
	for (bath, label) in ((discretefermionicbath([1.0], [0.5]; β=1.0), "μ=0"),
	                      (fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0, μ=0.7), "μ=0.7"))
		lat = GrassmannLattice(N=4, δτ=0.25, contour=:imag)
		kernel = correlationfunction(bath, lat)
		full = influence_functional(lat, kernel; alg=PartialIF()).state
		rowchain = GTEMPO.partial_if(lat, kernel; band=1)
		@test gd(full, rowchain) < 1e-12          # fullfactor vs 逐行因子积链
		exact = influence_functional(lat, kernel; alg=ExactTTI()).state
		imagchain = vacuumstate(lat)
		oldfac = Any[]
		GTEMPO.imag_factors(f -> push!(oldfac, f), spec, lat, kernel, 1, ExactTTI(), Dict{Symbol,Any}())
		for f in oldfac
			imagchain, _ = GM.mult(imagchain, f; trunc=Truncation(D=512, rtol=1e-14))
		end
		@test gd(exact, rowchain) < 1e-12         # ExactTTI(新) vs PartialIF（互换①）
		@test gd(exact, imagchain) < 1e-12        # tti_imag vs 旧 memory_cell 链
	end

	# ---- ① BCS fullfactor vs 旧 bcs_partial_if ----
	let latb = GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
		bathb = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0), Δ=0.4)
		kernelb = correlationfunction(bathb, latb)
		inew = influence_functional(latb, kernelb; alg=PartialIF()).state
		iold = GTEMPO.bcs_partial_if(latb, kernelb)
		@test gd(inew, iold) < 1e-12              # 反常象限整形式 vs 逐行链
	end

	# ---- ② 实时互换（连续谱阻尼核；偏差 = 拟合残差档，旧路径同源同 rtol）----
	let lat = GrassmannLattice(N=4, δt=0.15, contour=:real)
		spec_c = spectrum(ϵ -> sqrt(1 - (ϵ / 10)^2) / π, lb=-10, ub=10)
		bath = fermionicbath(spec_c, β=1.0)
		kernel = correlationfunction(bath, lat)
		exact = influence_functional(lat, kernel; alg=ExactTTI()).state
		rowchain = influence_functional(lat, kernel; alg=PartialIF()).state
		@test gd(exact, rowchain) < 1e-9          # cross_rtol=1e-10 拟合残差档
	end

	# ---- ② 实时增长极点守卫（DiracDelta |q|=1 边界极点响抛，现行语义；
	#       折叠朝向随链长/位型抖动轮盘——N=3 为本机实测 :right 位型）----
	let lat = GrassmannLattice(N=3, δt=0.15, contour=:real)
		bath = fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0)
		kernel = correlationfunction(bath, lat)
		@test_throws ArgumentError influence_functional(lat, kernel; alg=ExactTTI())
	end

	# ---- ② 混合互换（Pairing 直读路线——机器精度）----
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

	# ---- BCS + ExactTTI 拒绝（D7）----
	let latb = GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag, symmetry=Z2Symmetry())
		bathb = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0), Δ=0.4)
		kernelb = correlationfunction(bathb, latb)
		@test_throws ArgumentError influence_functional(latb, kernelb; alg=ExactTTI())
	end
end
