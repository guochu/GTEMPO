# BCS（D7 裁决：仅 :z2 对称性；反常通道 partial IF vs bcs_ed ED 参照）
println("========== Part 3: BCS (Z2 only) ==========")
@testset "BCS: Δ=0 退化为两条独立正常带" begin
	U = 0.0; ϵ_d = 0.8; ω = 1.0; α = 0.5
	β = 1.0; Nτ = 16; δτ = β / Nτ
	lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag, symmetry=Z2Symmetry())
	model = AndersonIM(U=U, ϵ_d=-ϵ_d)
	bathn = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	bath0 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=0)
	Ibcs = bcs_hybriddynamics(lat, bath0)
	kernel = correlationfunction(bathn, lat)
	In, _ = hybriddynamics(lat, kernel, PartialIF())
	K = sysdynamics(lat, model)
	Zb = integrate(K, Ibcs)
	Zn = integrate(K, In)
	@test value(integrate_scaled(K, Ibcs) / integrate_scaled(K, In)) ≈ 1.0 atol = 1e-8
end

@testset "BCS: 反常 IF vs bcs_ed（虚时 Z/occupation/GF）" begin
	U = 1.0; ϵ_ed = 0.8; ω = 1.0; α = 0.5
	β = 0.8; Nτ = 16; δτ = β / Nτ
	for Δ in (0.6, 0.2 + 0.45im)
		# ⚠️ ED 参照必须用 bcs_ed 自带算符（Ab/Bb）——其本地约定 occupied=首基、
		# 自带 JW 串，与 fermion_annihilators 的位基约定不同构（2026-11-16 教训：
		# 错配曾制造「occ 0.55 vs 0.32」假 bug——port 数值一直正确，反常通道对
		# Z/occ/奇 matsubara 全部 ED 吻合）
		H, Ab, Bb, H0 = bcs_ed(U, ϵ_ed; ω₀=ω, α=α, Δ=Δ)

		lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag, symmetry=Z2Symmetry())
		model = AndersonIM(U=U, ϵ_d=-ϵ_ed)
		bath = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
		Is = bcs_hybriddynamics(lat, bath)
		K = sysdynamics(lat, model)
		Z = integrate(K, Is)
		# BCS 浴自身带 Δ 配对项（破坏粒子数）→ 归一化分母 = 配对浴配分（复刻 bcs_ed 浴块）
		p1 = spin_half_matrices()
		n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
		I2 = one(n̂)
		Hbath = ω * (kron(n̂, I2) + kron(I2, n̂)) -
		        (Δ * kron(JW * σ₊, σ₊) - conj(Δ) * kron(JW * σ₋, σ₋))
		bathZ = tr(exp(-β * Hbath))
		@test real(Z) ≈ real(tr(exp(-β * H))) / bathZ rtol = 2e-2
		occ = occupation(lat, model, nothing, Is; j=1, branch=:τ, band=1, Z=Z)
		n_up = real(tr(Ab' * Ab * exp(-β * H))) / real(tr(exp(-β * H)))
		# 容差 = 虚时行因子化一阶 O(δτ) 常数（δτ = β/16）
		@test occ ≈ n_up atol = 1e-3
		# 奇插入穿反常 IF 的 matsubara（τ=β 点 = −⟨n⟩ 恒等式）
		m = matsubara(lat, model, nothing, Is; i=Nτ, j=0, Z=Z)
		@test m ≈ -n_up atol = 1e-3
	end
end

@testset "BCS: coupling_norm 语义（裸谱强度，Δ/μ/β 不进）" begin
	# Weyl 标量的 BCS 定义 = 每带裸谱强度（bcsweyl/README.md 推导：反常通道
	# 零贡献——Bogoliubov 正则性 u²+v²=1；接触缺口提取 cn*(Δ)/Σα ≡ 1 ≤ 4e-4）
	α = 0.5
	@test coupling_norm(bcsbath(fermionicbath(DiracDelta(ω=1.0, α=α), β=1.0), Δ=0.6)) == α
	@test coupling_norm(bcsbath(fermionicbath(DiracDelta(ω=1.0, α=α), β=1.0), Δ=0)) ==
	      coupling_norm(fermionicbath(DiracDelta(ω=1.0, α=α), β=1.0))   # Δ=0 与 normal 一致
	@test coupling_norm(bcsbath(fermionicbath(DiracDelta(ω=1.0, α=α), β=2.0, μ=0.3), Δ=1.0)) == α   # β/μ 不进
	vac = bcsvacuum(DiracDelta(ω=1.0, α=α), Δ=0.4)
	@test coupling_norm(vac) == α
	@test coupling_norm(discretebcsbath([1.0, 2.0], [α, α / 2]; β=1.0, Δ=0.3)) ≈ 1.5α
	# 带对 K 侧总量 = 2×（两带各 α；contact_extract 实证 cn* = 2α）
	bath = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=α), β=1.0), Δ=0.5)
	@test 2coupling_norm(bath) ≈ 1.0
end

@testset "BCS: 对称性守卫（D7 仅 :z2）" begin
	for sym in (U1Symmetry(), U1U1Symmetry(), Z2Z2Symmetry())
		@test_throws ArgumentError GTEMPO.bcs_symmetry_guard(sym)
	end
	lat = GrassmannLattice(N=4, δτ=0.25, bands=2, contour=:imag, symmetry=U1Symmetry())
	bath = bcsbath(fermionicbath(DiracDelta(ω=1.0, α=0.5), β=1.0), Δ=0.6)
	@test_throws ArgumentError bcs_hybriddynamics(lat, bath)
end
