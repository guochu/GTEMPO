@testset "API: integrate / integrate_scaled / ScaledScalar" begin
	# 多链 Berezin 迹：非交换链序只允许 (x, Is...) 形制；标度比值工具
	lat = GrassmannLattice(N=2, δτ=0.1, bands=1, contour=:imag)
	spec = GTEMPO.sectortype(Z2Symmetry())
	Is, _ = hybriddynamics(lat, correlationfunction(discretefermionicbath([1.0], [0.5]; β=lat.β), lat), PartialIF())
	K = sysdynamics(lat, ToulouseIM(ϵ_d=0.3))
	z = integrate_scaled(K, Is)
	@test z isa ScaledScalar
	@test value(z) ≈ integrate(K, Is) atol = 1.0e-12
	@test value(z / z) ≈ 1.0 atol = 1.0e-12
	# integrate(c::EnvironmentCache) 同值
	cache = environments_cache(K, Is)
	@test value(GM.Zvalue(cache)) ≈ integrate(K, Is) atol = 1.0e-12
end

@testset "API: mult / truncate! / truncerr 标准通道" begin
	lat = GrassmannLattice(N=2, δτ=0.1, bands=1, contour=:imag)
	spec = GTEMPO.sectortype(Z2Symmetry())
	spaces = site_spaces(spec, length(lat) ÷ 2, 1)
	x = randomgmps(ComplexF64, spaces; D=4)
	y = randomgmps(ComplexF64, spaces; D=4)
	# 朴素楔积（无截断）与 mult(trunc=Truncation(...)) 一致到截断水平
	exact = x * y
	out, err = mult(x, y; trunc=Truncation(D=typemax(Int), rtol=1e-14))
	@test length(out) == length(x)
	@test value(integrate_scaled(out) / integrate_scaled(exact)) ≈ 1.0 atol = 1.0e-10
	@test err isa Float64 && err >= 0
	# truncate! 原地 + 精确（大 D）→ 值不变；小 D → 键维受限
	c, err = truncate!(copy(exact); trunc=Truncation(D=typemax(Int), rtol=1e-14))
	@test value(integrate_scaled(c) / integrate_scaled(exact)) ≈ 1.0 atol = 1.0e-10
	_, err2 = truncate!(c; trunc=Truncation(D=2, rtol=1e-14))
	@test all(≤(2), bonddims(c))
	@test err2 isa Float64 && err2 >= 0
	# 零链回归（R2）：coeff=0 单项式（system 零核路径）经 truncate!/mult
	# 不破形——迹归零、键维不塌 0（新策略全零谱每扇区保 1 条）
	charges = cell_charges(spec, 1)
	z0 = monomial(spec, repeat(charges, 2); coeff=0)
	zt, _ = truncate!(z0; trunc=Truncation(D=8, rtol=1e-12))
	@test value(integrate_scaled(zt)) == 0.0
	zz, _ = mult(zt, x; trunc=Truncation(D=8, rtol=1e-12))
	@test value(integrate_scaled(zz)) == 0.0
	@test all(>(0), bonddims(zz))
end
