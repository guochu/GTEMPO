# 多带（≤2 带；Kanamori norb=1 = 2 带）ED 参照
println("========== Part 4: multi-orbital (≤2 bands) ==========")
@testset "Kanamori norb=1（2 带）：虚时 Z / occupation" begin
	μ = 0.3; U = 2.0; J = 0.4
	ω = 1.0; α = 0.5
	δτ = 0.1; Nτ = 8; β = Nτ * δτ
	H, _, _, H0 = kanamori_ed(U=U, J=J, norb=1, μ=μ, ω=ω, α=α)
	cs = fermion_annihilators(4)
	d1, d2 = cs[1], cs[2]

	lat = GrassmannLattice(N=Nτ, δτ=δτ, bands=2, contour=:imag, symmetry=:z2)
	bath = discretefermionicbath([ω], [α]; β=β)
	kernel = hybridization(lat, bath)
	model = KanamoriIM(U=U, J=J, norb=1, μ=μ)
	Is = ntuple(b -> hybriddynamics(lat, bath, PartialIF(); band=b), 2)
	K = sysdynamics(lat, model)
	Z = integrate(K, Is...)
	# FIXME(Kanamori Z 归一化, 2026-11-15)：J≠0 时 port Z/bathZ = 1.43 vs ED 3.08
	#（54%，δτ 减半不收敛）；AndersonIM U≠0 同路径仅 6% 一阶常数。occ/nn 比值观
	# 测机器精度 ✓（物理正确），绝对常数疑在 J 四算符项的槽分解标定（T1–T7 族）。
	# port 虚时 Z = 浴归一化值（对 ED 需除纯浴配分；O(δτ) 行因子化常数在内）
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
