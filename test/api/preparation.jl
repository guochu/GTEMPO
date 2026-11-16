# 制备对象（SystemFactorPreparation）两层缓存验收（plan
# measurement_cache_refactor §5 Phase B）：逐位性（验收 3）、键欠定义回归
# （验收 4）、扫描流 prep==fresh、守卫。双栈对拍（验收 2）由 observables.jl
# 既有 golden 锚维持（preparation=nothing = 现状路径，零行为变化）。
using LinearAlgebra: Diagonal
import GrassmannMPS as GM

@testset "API: OperatorInsertion 结构化 ==/hash（F4）" begin
	spec = GTEMPO.sectortype(U1Symmetry())
	qc = var_charge(spec, 1, 1)
	d = GTEMPO._fock_ops(1)[1]
	a = OperatorInsertion(:+, 3, d; charge=qc)
	@test a == OperatorInsertion(:+, 3, d; charge=qc)          # 新构造对象逐字段相等
	@test hash(a) == hash(OperatorInsertion(:+, 3, d; charge=qc))
	@test a != OperatorInsertion(:-, 3, d; charge=qc)          # branch
	@test a != OperatorInsertion(:+, 4, d; charge=qc)          # boundary
	@test a != OperatorInsertion(:+, 3, 2d; charge=qc)         # matrix
	@test a != OperatorInsertion(:+, 3, d)                     # charge（nothing）
	@test a != OperatorInsertion(:+, 3, d; charge=qc, side=:right)   # side
	@test a != OperatorInsertion(:+, 3, d; charge=GTEMPO.qneg(qc))
	# 同切位双插入键：元组 hash/== 走结构化定义
	k1 = (:+, 3, (a, OperatorInsertion(:+, 0, d'; charge=GTEMPO.qneg(qc))))
	k2 = (:+, 3, (OperatorInsertion(:+, 3, d; charge=qc), OperatorInsertion(:+, 0, d'; charge=GTEMPO.qneg(qc))))
	@test k1 == k2 && hash(k1) == hash(k2)
	@test Dict(k1 => 1)[k2] == 1
end

@testset "API: 制备对象逐位性（验收 3：缓存命中 = 重算）" begin
	μ = 0.4; β = 0.8
	lat = GrassmannLattice(N=3, δt=0.1, bands=1, contour=:real, symmetry=U1Symmetry())
	model = ToulouseIM(ϵ_d=μ)
	q = exp(-β * μ)
	ρ0 = Matrix{ComplexF64}(Diagonal([1.0, q])); ρ0 ./= tr(ρ0)
	spec = GTEMPO.sectortype(U1Symmetry())
	qc = var_charge(spec, 1, 1)
	d = GTEMPO._fock_ops(1)[1]; ddag = d'
	greater_ins(i) = [OperatorInsertion(:+, i, d; charge=GTEMPO.qneg(qc)),
	                  OperatorInsertion(:+, 0, ddag; charge=qc)]
	prep = prepare_system(lat, model; order=1, initial_state=ρ0)
	# 无插入 + 奇插入对（扫描首/中/末点）：prep vs 无缓存重建逐位
	@test GM.todense(sysdynamics(prep)) == GM.todense(sysdynamics(lat, model; initial_state=ρ0))
	for i in (1, 2, 3)
		Kp = sysdynamics(prep; insertions=greater_ins(i))
		Kf = sysdynamics(lat, model; initial_state=ρ0, insertions=greater_ins(i))
		@test GM.todense(Kp) == GM.todense(Kf)
		# 命中路径：重复调用逐位不变（别名污染回归：缓存列不得被 append! 改写）
		@test GM.todense(sysdynamics(prep; insertions=greater_ins(i))) == GM.todense(Kp)
	end
	# 虚时轮廓（:τ 传播子键）同款逐位
	lati = GrassmannLattice(N=3, δτ=0.25, contour=:imag, symmetry=U1Symmetry())
	prepi = prepare_system(lati, model)
	Kip = sysdynamics(prepi; insertions=[OperatorInsertion(:τ, 2, d; charge=GTEMPO.qneg(qc)),
	                                     OperatorInsertion(:τ, 0, ddag; charge=qc)])
	Kif = sysdynamics(lati, model; insertions=[OperatorInsertion(:τ, 2, d; charge=GTEMPO.qneg(qc)),
	                                           OperatorInsertion(:τ, 0, ddag; charge=qc)])
	@test GM.todense(Kip) == GM.todense(Kif)
end

@testset "API: 缓存键完备性（验收 4：同格换插入产出不同且各自正确）" begin
	μ = 0.4; β = 0.8
	lat = GrassmannLattice(N=3, δt=0.1, bands=1, contour=:real, symmetry=U1Symmetry())
	model = ToulouseIM(ϵ_d=μ)
	q = exp(-β * μ)
	ρ0 = Matrix{ComplexF64}(Diagonal([1.0, q])); ρ0 ./= tr(ρ0)
	spec = GTEMPO.sectortype(U1Symmetry())
	qc = var_charge(spec, 1, 1)
	d = GTEMPO._fock_ops(1)[1]; ddag = d'
	prep = prepare_system(lat, model; order=1, initial_state=ρ0)
	n = GTEMPO._number_op(1, 1)
	# 同格 (:+,3)：n̂ / 奇对 / :right 侧 / 空切位——各自与无缓存重建逐位一致
	ins_n = [OperatorInsertion(:+, 3, n)]
	ins_dn = [OperatorInsertion(:+, 3, d; charge=GTEMPO.qneg(qc)),
	          OperatorInsertion(:+, 0, ddag; charge=qc)]
	ins_r = [OperatorInsertion(:+, 3, n; side=:right)]
	for ins in (ins_n, ins_dn, ins_r)
		@test GM.todense(sysdynamics(prep; insertions=ins)) ==
		      GM.todense(sysdynamics(lat, model; initial_state=ρ0, insertions=ins))
	end
	# 且互不相同（键漏字段 = 同值假命中）
	dn = GM.todense(sysdynamics(prep; insertions=ins_dn))
	nn_ = GM.todense(sysdynamics(prep; insertions=ins_n))
	rr = GM.todense(sysdynamics(prep; insertions=ins_r))
	ee = GM.todense(sysdynamics(prep))
	@test dn != nn_ && dn != ee && nn_ != ee && rr != nn_
end

@testset "API: 扫描流 prep==fresh（greater/lesser/occupation/current）" begin
	μ = 0.7; β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)
	lat = GrassmannLattice(N=4, δt=0.15, bands=2, contour=:real, symmetry=U1Symmetry())
	q = exp(-β * μ)
	ρ0 = Matrix{ComplexF64}(Diagonal([1.0, q, q, q * q])); ρ0 ./= tr(ρ0)
	model = AndersonIM(U=0.0, ϵ_d=μ)
	kernel = hybridization(lat, bath; coupling_norm=1.0)
	Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b, coupling_norm=1.0), 2)
	K = sysdynamics_weyl2order(lat, model; coupling_norm=1.0, initial_state=ρ0)
	zs = integrate_scaled(K, Is...)
	Z = exp(zs.logabs) * zs.phase
	prep = prepare_system(lat, model; order=2, coupling_norm=1.0, initial_state=ρ0)
	# greater/lesser 全点：prep 与逐点显式 kwargs（现状路径）标量逐位相等
	for i in 1:4
		@test greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, preparation=prep) ==
		      greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2, coupling_norm=1.0, initial_state=ρ0)
		@test lesser(lat, model, kernel, Is...; i=i, j=0, Z=Z, preparation=prep) ==
		      lesser(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2, coupling_norm=1.0, initial_state=ρ0)
	end
	@test occupation(lat, model, kernel, Is...; j=2, branch=:+, Z=Z, preparation=prep) ==
	      occupation(lat, model, kernel, Is...; j=2, branch=:+, Z=Z, order=2, coupling_norm=1.0, initial_state=ρ0)
	# 电流（k 固定 j 扫描形态，distinct 格最多）：prep == fresh
	@test electriccurrent(lat, model, kernel, Is...; k=3, band=1, Z=Z, preparation=prep) ==
	      electriccurrent(lat, model, kernel, Is...; k=3, band=1, Z=Z, order=2, coupling_norm=1.0, initial_state=ρ0)
	# clear_cache! 后重算仍一致
	clear_cache!(prep)
	@test greater(lat, model, kernel, Is...; i=2, j=0, Z=Z, preparation=prep) ==
	      greater(lat, model, kernel, Is...; i=2, j=0, Z=Z, order=2, coupling_norm=1.0, initial_state=ρ0)
	# lattice 不匹配守卫（δt 不同 → 传播子键失效域）
	lat2 = GrassmannLattice(N=4, δt=0.1, bands=2, contour=:real, symmetry=U1Symmetry())
	@test_throws ArgumentError greater(lat2, model, kernel, Is...; i=1, Z=Z, preparation=prep)
	# weyl2order(prep)：order 守卫 + 制备绑定 coupling_norm（逐位段用小链，
	# todense 合规）
	@test_throws ArgumentError sysdynamics_weyl2order(prepare_system(lat, model))
	lats = GrassmannLattice(N=3, δt=0.1, bands=1, contour=:real, symmetry=U1Symmetry())
	ρ0s = Matrix{ComplexF64}(Diagonal([1.0, exp(-β * μ)])); ρ0s ./= tr(ρ0s)
	prep1 = prepare_system(lats, ToulouseIM(ϵ_d=μ); order=2, coupling_norm=0.5, initial_state=ρ0s)
	@test GM.todense(sysdynamics_weyl2order(prep1)) ==
	      GM.todense(sysdynamics_weyl2order(lats, ToulouseIM(ϵ_d=μ); coupling_norm=0.5, initial_state=ρ0s))
end
