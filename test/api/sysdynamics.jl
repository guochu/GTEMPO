using LinearAlgebra: Diagonal, norm

@testset "API: sysdynamics（K 构造）" begin
	# 虚时 order=1：链长/站数/物理性
	lat = GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag)
	model = AndersonIM(U=1.0, ϵ_d=-0.5)
	K = sysdynamics(lat, model)
	@test K isa GM.GMPS && length(K) == length(lat) ÷ 2   # pair 布局：每胞 bands 站
	@test iscanonical(K)
	@test_throws UndefKeywordError sysdynamics_weyl2order(lat, model)   # 缺 coupling_norm kwarg
	@test sysdynamics(lat, model; order=2) isa GM.GMPS   # cn=0 合法（仅无 Weyl 标量）

	# 实时 order=1 / initial_state（实时专属）
	latr = GrassmannLattice(N=4, δt=0.05, bands=2, contour=:real)
	q = exp(-1.0 * (-0.5))
	ρ0 = Matrix{ComplexF64}(Diagonal([1.0, q, q, q * q] ./ (1 + 2q + q * q)))
	Kr = sysdynamics(latr, model; initial_state=ρ0)
	@test length(Kr) == length(latr) ÷ 2 && norm(Kr) > 0

	# 偶算符插入（charge=nothing）合法
	spec = symmetry_spec(:z2)
	d = GTEMPO._fock_ops(2)[1]
	ins_even = [OperatorInsertion(:+, 1, d)]
	K2 = sysdynamics(latr, model; insertions=ins_even, initial_state=ρ0)
	@test length(K2) == length(latr) ÷ 2
	# 单奇插入净荷非零 → 显式拒绝；奇配对（ddag@0 + d@1）合法
	latu1 = GrassmannLattice(N=4, δt=0.05, bands=2, contour=:real, symmetry=:u1)
	qc = GTEMPO.var_charge(symmetry_spec(:u1), 1, 2)
	ins_odd = [OperatorInsertion(:+, 1, d; charge=GTEMPO.qneg(qc))]
	@test_throws ArgumentError sysdynamics(latu1, model; insertions=ins_odd, initial_state=ρ0)
	ddag = d'
	ins_pair = [OperatorInsertion(:+, 0, ddag; charge=qc),
	            OperatorInsertion(:+, 1, d; charge=GTEMPO.qneg(qc))]
	Ku1 = sysdynamics(latu1, model; insertions=ins_pair, initial_state=ρ0)
	@test length(Ku1) == length(latu1) ÷ 2
	@test Ku1 isa GM.GMPS

	# 零时轮廓 / 单格链
	lat1 = GrassmannLattice(N=1, δτ=0.1, contour=:imag)
	K1 = sysdynamics(lat1, ToulouseIM(ϵ_d=-0.5))
	@test length(K1) == 1
end

@testset "API: sysdynamics_weyl2order" begin
	lat = GrassmannLattice(N=3, δt=0.1, bands=2, contour=:real)
	model = AndersonIM(U=0.5, ϵ_d=-0.5)
	K = sysdynamics_weyl2order(lat, model; coupling_norm=1.0)
	@test length(K) == length(lat) ÷ 2 && norm(K) > 0
	# Weyl 标量 = exp(cn/4·Σ step²)：order=1 链 × e^{−Nt·dt²/4·...} 平移不改变迹，
	# 这里只验证量级与键维有限
	@test all(>(0), bond_dimensions(K))
	# K(weyl, order=2) ∧ 自身 = K₂₂ 因子化检验（U=0 精确路径在 ED 电池覆盖）
end
