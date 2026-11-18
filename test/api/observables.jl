# 观测量电池：inserted-K 精确路线（golden 对拍 U1GTEMPO2 机器精度锚）+
# ED 参照（单模浴解析动力学）+ 算符路线（cached_*）+ 电流。
using LinearAlgebra: Diagonal, eigen

@testset "golden: greater/lesser vs U1GTEMPO2（机器精度锚）" begin
	μ = 0.7; β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)
	lat = GrassmannLattice(N=4, δt=0.15, bands=2, contour=:real, symmetry=U1Symmetry())
	q = exp(-β * μ)
	ρ0 = Diagonal(ComplexF64[1.0, q, q, q * q]); ρ0 /= tr(ρ0)
	model = AndersonIM(U=0.0, ϵ_d=μ)
	kernel = hybridization(lat, bath; coupling_norm=1.0)
	Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b, coupling_norm=1.0), 2)
	K = sysdynamics_weyl2order(lat, model; coupling_norm=1.0, initial_state=ρ0)
	zs = integrate_scaled(K, Is...)
	Z = exp(zs.logabs) * zs.phase
	@test Z ≈ 0.9998529294302793 atol = 1e-10
	# 参照栈 solve(greater) 全序列（U1GTEMPO2;机器精度锚）
	ref_gt = [-0.06945314205699218 - 0.6608068990426897im,
	          -0.13549191017047374 - 0.6390680792460552im,
	          -0.1948904437761824 - 0.6041604002540936im,
	          -0.2447894924470817 - 0.5579925365001444im]
	ref_lt = [0.03448940964763923 + 0.328146994708573im,
	          0.06728329137328251 + 0.3173518162757938im,
	          0.09677973022863783 + 0.3000171758989178im,
	          0.12155886447177881 + 0.27709089325126396im]
	for i in 1:4
		g = greater(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2, coupling_norm=1.0, initial_state=ρ0)
		l = lesser(lat, model, kernel, Is...; i=i, j=0, Z=Z, order=2, coupling_norm=1.0, initial_state=ρ0)
		@test g ≈ ref_gt[i] atol = 5e-14
		@test l ≈ ref_lt[i] atol = 5e-14
	end
	# occupation = 总占据（两带对称 ⟹ 2× 单带值）
	occ = occupation(lat, model, kernel, Is...; j=1, branch=:+, Z=Z, order=2,
	                 coupling_norm=1.0, initial_state=ρ0)
	@test occ ≈ 2 * 0.33111020348052433 atol = 5e-14
end

@testset "golden: matsubara vs U1GTEMPO2（虚时）" begin
	μ = 0.7; β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)
	lat = GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag, symmetry=U1Symmetry())
	model = AndersonIM(U=0.0, ϵ_d=μ)
	kernel = hybridization(lat, bath)
	Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b), 2)
	K = sysdynamics(lat, model)
	Z = integrate(K, Is...)
	@test Z ≈ 2.543041612650098 atol = 1e-10
	ref_m = [-0.532286550768878, -0.44990671645915964, -0.3916734863440781, -0.3499217507041494]
	for i in 1:4
		m = matsubara(lat, model, kernel, Is...; i=i, j=0, Z=Z)
		@test m ≈ ref_m[i] atol = 5e-14
	end
	# 平移不变（平衡态）：⟨d(τ4−τ3)d†(0)⟩ = ⟨d(τ1−0)d†(0)⟩ ⟹ m43 = ref_m[1]
	m43 = matsubara(lat, model, kernel, Is...; i=4, j=3, Z=Z)
	@test m43 ≈ ref_m[1] atol = 1e-12
	m34 = matsubara(lat, model, kernel, Is...; i=3, j=4, Z=Z)
	@test isfinite(real(m34))
end

@testset "ED: 单模浴 occupation / greater / lesser / matsubara" begin
	μ = 0.4; β = 0.8
	ω, α = 1.0, 0.5
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0.0, bathspecs=[(ω, α)], nbands=1)

	# 虚时
	lat = GrassmannLattice(N=8, δτ=β/8, contour=:imag)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = hybridization(lat, bath)
	Is = hybriddynamics(lat, bath, PartialIF(); band=1)
	K = sysdynamics(lat, ToulouseIM(ϵ_d=μ))
	Z = integrate(K, Is)
	# Grassmann 路线 Z = 浴归一化（缺纯浴配分因子，Z_B 铁律）+ 行因子化 O(δτ) 常数
	@test Z ≈ tr(exp(-β * H)) / (1 + exp(-β * ω)) atol = 3e-2
	τs = [i * lat.δτ for i in 1:lat.k]
	gτ_ref = correlation_2op_1τ(H, a, adag, τs; β=β)
	for i in 1:lat.k
		m = matsubara(lat, ToulouseIM(ϵ_d=μ), kernel, Is; i=i, j=0, Z=Z)
		@test m ≈ -gτ_ref[i] atol = 1e-3      # 时序因子 (i ≥ 0) = −⟨a(τ)a†⟩
	end

	# 实时（热初态 = 杂质热态 × 浴真空；H0 = 退耦 Hamiltonian）
	lat = GrassmannLattice(N=3, δt=0.05, contour=:real)
	kernel = hybridization(lat, bath)
	Is = hybriddynamics(lat, bath, PartialIF(); band=1)
	d_imp = fermion_annihilators(1)[1]
	himp = μ * d_imp' * d_imp
	ρ0 = Matrix{ComplexF64}(exp(-β * himp)); ρ0 /= tr(ρ0)
	model = ToulouseIM(ϵ_d=μ)
	K = sysdynamics(lat, model; initial_state=ρ0)
	Z = integrate(K, Is)
	ts = [i * lat.δt for i in 1:lat.k]
	gt_ref, lt_ref = greater_lesser_ed(H, a, adag, H0, ts, β)
	for i in 1:lat.k
		g = greater(lat, model, kernel, Is; i=i, j=0, Z=Z, initial_state=ρ0)
		l = lesser(lat, model, kernel, Is; i=i, j=0, Z=Z, initial_state=ρ0)
		@test g ≈ gt_ref[i] atol = 1e-5
		@test l ≈ lt_ref[i] atol = 1e-5
	end

	# U≠0（两带 U·n₁n₂ + 每带耦合两条浴模——与 port 核的 Σ_m 求和约定一致）
	H2 = zeros(Float64, 16, 16)
	ed_density!(H2, 1, μ); ed_density!(H2, 2, μ)
	ed_density_interaction!(H2, 1, 2, 1.0)
	H02 = copy(H2)
	for m in 1:2
		ed_density!(H2, 2 + m, ω); ed_density!(H02, 2 + m, ω)
	end
	for b in 1:2, m in 1:2
		v = sqrt(α)
		ed_tunneling!(H2, 2 + m, b, v); ed_tunneling!(H2, b, 2 + m, v)
	end
	a2 = fermion_annihilators(4)[1]
	lat2 = GrassmannLattice(N=6, δτ=β/6, bands=2, contour=:imag)
	bath2 = discretefermionicbath([ω, ω], [α, α]; β=β)
	kernel2 = hybridization(lat2, bath2)
	Is2 = ntuple(b -> hybriddynamics(lat2, bath2, PartialIF(); band=b), 2)
	K2 = sysdynamics(lat2, AndersonIM(U=1.0, ϵ_d=μ))
	Z2 = integrate(K2, Is2...)
	bathZ2 = (1 + exp(-β * ω))^2
	# 行因子化 IF 的一阶 O(δτ) 常数（实测 ~0.45/δτ，δτ 减半误差减半；golden 锚为
	# port≡参照离散化的机器精度，ED 比较的容差按此常数标定）
	@test real(Z2) ≈ tr(exp(-β * H2)) / bathZ2 rtol = 8e-2
	occ = occupation(lat2, AndersonIM(U=1.0, ϵ_d=μ), kernel2, Is2...; j=1, branch=:τ, band=1, Z=Z2)
	n_ref = real(tr(a2' * a2 * exp(-β * H2))) / real(tr(exp(-β * H2)))
	@test occ ≈ n_ref atol = 5e-3
end

@testset "ED: nn 密度-密度（时序关联）" begin
	μ = 0.4; β = 0.8; ω, α = 1.0, 0.5
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0.0, bathspecs=[(ω, α)], nbands=1)
	lat = GrassmannLattice(N=3, δt=0.05, contour=:real)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = hybridization(lat, bath)
	Is = hybriddynamics(lat, bath, PartialIF(); band=1)
	d_imp = fermion_annihilators(1)[1]
	himp = μ * d_imp' * d_imp
	ρ0 = Matrix{ComplexF64}(exp(-β * himp)); ρ0 /= tr(ρ0)
	model = ToulouseIM(ϵ_d=μ)
	K = sysdynamics(lat, model; initial_state=ρ0)
	Z = integrate(K, Is)
	n = adag * a
	ts = [i * lat.δt for i in 1:lat.k]
	nn_ref = correlation_2op_1t(H, n, n, Matrix{ComplexF64}(exp(-β * H0)), ts)
	for i in 1:lat.k
		v = nn(lat, model, kernel, Is; i=i, j=1, b1=:+, b2=:+, Z=Z, initial_state=ρ0)
		# 边界插入点接触项 O(δt)（与 greater 同源；旧栈同款量级）
		@test v ≈ nn_ref[i] atol = 3e-3
	end
	# 同格点退化：nn(i,i) = occupation(i)
	o = occupation(lat, model, kernel, Is; j=2, branch=:+, Z=Z, initial_state=ρ0)
	@test nn(lat, model, kernel, Is; i=2, j=2, b1=:+, b2=:+, Z=Z, initial_state=ρ0) ≈ o atol = 1e-12
end

@testset "算符路线：cached_* vs inserted-K（D5 双轨）" begin
	μ = 0.4; β = 0.8; ω, α = 1.0, 0.5
	lat = GrassmannLattice(N=3, δt=0.05, contour=:real)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = hybridization(lat, bath)
	Is = hybriddynamics(lat, bath, PartialIF(); band=1)
	d_imp = fermion_annihilators(1)[1]
	himp = μ * d_imp' * d_imp
	ρ0 = Matrix{ComplexF64}(exp(-β * himp)); ρ0 /= tr(ρ0)
	model = ToulouseIM(ϵ_d=μ)
	K = sysdynamics(lat, model; initial_state=ρ0)
	Z = integrate(K, Is)
	cache = environments_cache(K, Is)
	@test cached_Z(cache) ≈ Z atol = 1e-10
	# 参照点 j=1（j=0 是 ρ₀ 边界矩阵切位，无变量位，单子式路线不可达）
	for i in 2:lat.k
		g_ins = greater(lat, model, kernel, Is; i=i, j=1, Z=Z, initial_state=ρ0)
		# 算符路线 = 胞坐标 0⁺ 点原始关联（无 −i 因子），与边界插入差 O(接触项)
		#（D5 裁决，容差按胞宽标定）
		g_cached = cached_correlation(lat, K, Is; i=i, j=1, kind=:greater)
		@test abs(g_cached - g_ins / (-1.0im)) < 5e-3
	end
	# lattice_monomial 荷闭合（SparseGMPS 支撑 = 触碰站闭包）
	idx = (GrassmannContourIndex(2, :+, 1, false), GrassmannContourIndex(1, :+, 1, true))
	op = lattice_monomial(lat, idx)
	@test op isa SparseGMPS
end

@testset "电流：raw_lead_correlator / electriccurrent / heatcurrent" begin
	μ = 0.7; β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)
	lat = GrassmannLattice(N=4, δt=0.15, bands=2, contour=:real, symmetry=U1Symmetry())
	q = exp(-β * μ)
	ρ0 = Diagonal(ComplexF64[1.0, q, q, q * q]); ρ0 /= tr(ρ0)
	model = AndersonIM(U=0.0, ϵ_d=μ)
	kernel = hybridization(lat, bath; coupling_norm=1.0)
	Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b, coupling_norm=1.0), 2)
	K = sysdynamics_weyl2order(lat, model; coupling_norm=1.0, initial_state=ρ0)
	Z = integrate(K, Is...)
	# raw_lead_correlator = lesser 图案（⟨d†(+,k) d(b2,j)⟩；同切位 a†a）
	r = raw_lead_correlator(lat, model, kernel, Is...; k=2, j=2, b2=:+, band=1, Z=Z,
	                        order=2, coupling_norm=1.0, initial_state=ρ0)
	@test r ≈ 0.329036909179571 atol = 1e-12          # = n(+,2)（ED 值锚）
	r2 = raw_lead_correlator(lat, model, kernel, Is...; k=2, j=1, b2=:-, band=1, Z=Z,
	                         order=2, coupling_norm=1.0, initial_state=ρ0)
	@test r2 ≈ 0.6621914237309042 + 0.06963043465201318im atol = 1e-12
	# 电流有限实值（U=0 平衡态净电流 ≈ 0；权重公式合法性）
	I = electriccurrent(lat, model, kernel, Is...; k=3, band=1, Z=Z,
	                    order=2, coupling_norm=1.0, initial_state=ρ0)
	@test isfinite(real(I)) && isfinite(imag(I))
	# 平衡态净电流 = 0，离散化伪影 O(δt)（旧栈 ED 电流对比同款量级）
	@test abs(I) < 5e-2
	# 热流 = ω 加权谱的电流；平衡态同理为零
	Q = heatcurrent(lat, model, bath, Is...; k=3, band=1, Z=Z,
	                order=2, coupling_norm=1.0, initial_state=ρ0)
	@test isfinite(real(Q)) && abs(Q) < 5e-2
end

@testset "电流快速版：流算符 GMPS（两轨，备份 electriccurrent_fast 同型）" begin
	μ = 0.7; β = 1.0
	bath = discretefermionicbath([1.0], [0.5]; β=β)
	lat = GrassmannLattice(N=4, δt=0.15, bands=2, contour=:real, symmetry=U1Symmetry())
	q = exp(-β * μ)
	ρ0 = Diagonal(ComplexF64[1.0, q, q, q * q]); ρ0 /= tr(ρ0)
	model = AndersonIM(U=0.0, ϵ_d=μ)
	kernel = hybridization(lat, bath; coupling_norm=1.0)
	Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b, coupling_norm=1.0), 2)
	K = sysdynamics_weyl2order(lat, model; coupling_norm=1.0, initial_state=ρ0)
	zs = integrate_scaled(K, Is...)
	Z = exp(zs.logabs) * zs.phase
	cache = environments_cache(K, Is...)
	ηpp = Matrix(kernel.blocks[(:+, :+)])
	ηpm = -Matrix(kernel.blocks[(:+, :-)])
	# ① 机器精度：两轨快速版（快照扫掠 / J∧K 单发积分）== 逐项 cached_gf
	#   线性和（cached_gf = 原始收缩，除 Z 后同口径——SparseGMPS 直和的线性）
	for k in 2:4, band in (1, 2)
		term = 0.0im
		row = GrassmannContourIndex(Int(k), :+, band, true)
		for j in 1:k-1, (b2, η) in ((:+, ηpp), (:-, ηpm))
			term += η[k, j] * cached_gf(lat, (row, GrassmannContourIndex(Int(j), b2, band, false)), K, Is...)
		end
		linear = 2 * term / (Z * lat.δt)
		fast_c = cached_electriccurrent_fast(lat, kernel, K, Is...; k=k, band=band, cache=cache)
		fast_d = electriccurrent_fast(lat, model, kernel, Is...; k=k, band=band, Z=Z,
		                              order=2, coupling_norm=1.0, initial_state=ρ0)
		@test fast_c ≈ linear atol = 1e-12
		@test fast_d ≈ linear atol = 1e-12
	end
	# ② k = 1（空和）短路零电流
	@test cached_electriccurrent_fast(lat, kernel, K, Is...; k=1) == 0.0im
	@test electriccurrent_fast(lat, model, kernel, Is...; k=1, Z=Z) == 0.0im
	# ③ D5 跨语义：快速版 = 胞坐标流算符（备份 fast 同款），与边界插入 slow 差
	#   O(胞内接触项)，经 2/δt·Σ|η| 放大 —— 实测 δt=0.15 时 |Δ| ≈ 0.054·k/2 档，
	#   容差按此标定（0.08k）；边界语义锚定路径仍是 electriccurrent
	for k in 2:4
		slow = electriccurrent(lat, model, kernel, Is...; k=k, band=1, Z=Z,
		                       order=2, coupling_norm=1.0, initial_state=ρ0)
		fast = cached_electriccurrent_fast(lat, kernel, K, Is...; k=k, band=1, cache=cache)
		@test abs(fast - slow) < 0.08 * k
	end
end
