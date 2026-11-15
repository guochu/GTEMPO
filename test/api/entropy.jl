@testset "API: entanglement entropy（GM 引擎面）" begin
	spec = symmetry_spec(:z2)
	spaces = site_spaces(spec, 8, 1)
	x = randomgmps(ComplexF64, spaces; D=8)
	canonicalize!(x)

	for bond in 1:7
		p = entanglement_spectrum(x; bond=bond)
		@test length(p) == bonddims(x)[bond]
		@test sum(p) ≈ 1
		@test all(p .>= -1e-14)
		@test schmidt_values(x; bond=bond) ≈ sqrt.(p)
		@test entanglement_entropy(x; bond=bond) ≈
		      -sum(q * log(q) for q in p if q > 0)
		@test entanglement_entropy(x; bond=bond, α=2) ≈ -log(sum(abs2, p))
	end
	@test entanglement_entropy(x) ≈ entanglement_entropy(x; bond=4)
	@test_throws BoundsError entanglement_entropy(x; bond=0)
	@test_throws BoundsError entanglement_entropy(x; bond=8)

	# Schmidt 值未初始化时就地 canonicalize
	y = deepcopy(x)
	GM.unset_svectors!(y)
	@test GM.svectors_uninitialized(y)
	S = entanglement_entropy(y; bond=4)
	@test !GM.svectors_uninitialized(y)
	@test S ≈ entanglement_entropy(x; bond=4)

	# 真空态（无纠缠）各键熵为零
	v = vacuumstate(GrassmannLattice(N=8, δτ=0.1, contour=:imag))
	for bond in 1:7
		@test entanglement_entropy(v; bond=bond) ≈ 0 atol=1e-14
		@test entanglement_spectrum(v; bond=bond) ≈ [1.0]
	end
end
