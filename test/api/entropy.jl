using LinearAlgebra: normalize!

@testset "API: entanglement entropy" begin
	L, D = 8, 8
	x = randomgmps(L, D=D)
	canonicalize!(x, alg=Orthogonalize(trunc=NoTruncation(), normalize=false))

	for bond in 1:L-1
		p = entanglement_spectrum(x; bond=bond)
		@test length(p) == bonddim(x, bond)
		@test sum(p) ≈ 1
		@test all(p .>= -1e-14)
		@test schmidt_values(x; bond=bond) ≈ sqrt.(p)
		@test entanglement_entropy(x; bond=bond) ≈
		      -sum(q * log(q) for q in p if q > 0)
		# 二阶 Rényi 熵：-log(tr ρ²) = -log(Σ p²)
		@test entanglement_entropy(x; bond=bond, α=2) ≈ -log(sum(abs2, p))
	end

	# 默认 bond = 链中央
	@test entanglement_entropy(x) ≈ entanglement_entropy(x; bond=div(L, 2))

	# bond 越界
	@test_throws BoundsError entanglement_entropy(x; bond=0)
	@test_throws BoundsError entanglement_entropy(x; bond=L)

	# Schmidt 值未初始化时就地 canonicalize（与 FMA 的语义一致）
	y = deepcopy(x)
	unset_svectors!(y)
	@test svectors_uninitialized(y)
	S = entanglement_entropy(y; bond=4)
	@test !svectors_uninitialized(y)
	@test S ≈ entanglement_entropy(x; bond=4)

	# 真空态（无纠缠）各键熵为零
	vac = GrassmannMPS(Float64, L)
	for bond in 1:L-1
		@test entanglement_entropy(vac; bond=bond) ≈ 0 atol=1e-14
		@test entanglement_spectrum(vac; bond=bond) ≈ [1.0]
	end
end
