# 多带（≤2 带；Kanamori norb=1 = 2 带）ED 参照
println("========== Part 4: multi-orbital (≤2 bands) ==========")
@testset "Kanamori norb=1（2 带）：虚时 Z / occupation" begin
	μ = 0.3; U = 2.0; J = 0.4
	ω = 1.0; α = 0.5
	δτ = 0.1; Nτ = 8; β = Nτ * δτ
	H, _, _, H0 = kanamori_ed(U=U, J=J, norb=1, μ=μ, ω=ω, α=α)
	cs = fermion_annihilators(4)
	d1, d2 = cs[1], cs[2]

	lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag, symmetry=Z2Symmetry())
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = correlationfunction(bath, lat)
	model = KanamoriIM(U=U, J=J, norb=1, μ=μ)
	Is, _ = hybriddynamics(lat, kernel, PartialIF())
	K = sysdynamics(lat, model)
	Z = integrate(K, Is...)
	# port 虚时 Z = 浴归一化配分（Z_full/Z_B——Z_B 已在归一化关联函数处除掉），
	# ED 对比口径 = ED_abs/Z_B（(1+e^{-βω})² = 双带纯浴配分）。实测比值
	# 0.974（Nτ=8）→ 0.987（Nτ=16）一阶收敛 ✓。
	# 勘误（2026-11-16 上游裁决，GrassmannMPS/bugreports/
	# 2026-11-16-upstream-verdict-multiband-z.md）：本文件 2026-11-15 的
	# 「J≠0 缺 54%、疑 J 四算符槽分解标定」记录 = 对 port 侧重复除 Z_B 的双除
	# 伪象（Z_port/Z_B = 1.43 误配 ED/Z_B = 3.08）——occ/nn 机器精度与 J 槽
	# 分解标定本就无缺陷，与撤回的 bands=4 绝对 Z 报告同一类口径错。
	@test real(Z) ≈ real(tr(exp(-β * H))) / (1 + exp(-β * ω))^2 rtol = 4e-2
	for band in 1:2
		occ = occupation(lat, model, kernel, Is...; j=1, branch=:τ, band=band, Z=Z)
		d = band == 1 ? d1 : d2
		n_ed = real(tr(d' * d * exp(-β * H))) / real(tr(exp(-β * H)))
		@test occ ≈ n_ed atol = 1e-3
	end
	# nn（band 内 密度-密度）
	nn1 = nn(lat, model, kernel, Is...; i=1, j=1, b1=:τ, b2=:τ, band=(1, 1), Z=Z)
	n1 = real(tr(d1' * d1 * exp(-β * H))) / real(tr(exp(-β * H)))
	@test nn1 ≈ n1 atol = 1e-3   # 同格同带 n² = n
end
