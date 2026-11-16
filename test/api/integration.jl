@testset "API: integrate / integrate_scaled / ScaledScalar" begin
	# 多链 Berezin 迹：非交换链序只允许 (x, Is...) 形制；标度比值工具
	lat = GrassmannLattice(N=2, δτ=0.1, bands=1, contour=:imag)
	spec = symmetry_spec(:z2)
	Is = hybriddynamics(lat, discretefermionicbath([1.0], [0.5]; β=lat.β), PartialIF(); band=1)
	K = sysdynamics(lat, ToulouseIM(ϵ_d=0.3))
	z = integrate_scaled(K, Is)
	@test z isa ScaledScalar
	@test unscale(z) ≈ integrate(K, Is) atol = 1.0e-12
	@test scaled_ratio(z, z) ≈ 1.0 atol = 1.0e-12
	# integrate(c::EnvironmentCache) 同值
	cache = environments_cache(K, Is)
	@test unscale(GM.Zvalue(cache)) ≈ integrate(K, Is) atol = 1.0e-12
end

@testset "API: multiply / compress! / 截断记账" begin
	lat = GrassmannLattice(N=2, δτ=0.1, bands=1, contour=:imag)
	spec = symmetry_spec(:z2)
	spaces = site_spaces(spec, length(lat) ÷ 2, 1)
	x = randomgmps(ComplexF64, spaces; D=4)
	y = randomgmps(ComplexF64, spaces; D=4)
	# 朴素楔积（无截断）与 multiply(trunc=Truncation(...)) 一致到截断水平
	exact = x * y
	reports = CompressionReport[]
	out = multiply(x, y; trunc=Truncation(maxdim=typemax(Int), rtol=1e-14), reports=reports)
	@test length(out) == length(x)
	@test scaled_ratio(integrate_scaled(out), integrate_scaled(exact)) ≈ 1.0 atol = 1.0e-10
	@test length(reports) == 1 && reports[1] isa CompressionReport
	# compress! 原地 + 精确（大 D）→ 值不变；小 D → 键维受限
	c = copy(exact)
	r = compress!(c; trunc=Truncation(maxdim=typemax(Int), rtol=1e-14))
	@test scaled_ratio(integrate_scaled(c), integrate_scaled(exact)) ≈ 1.0 atol = 1.0e-10
	r2 = compress!(c; trunc=Truncation(maxdim=2, rtol=1e-14))
	@test all(≤(2), bond_dimensions(c))
end
