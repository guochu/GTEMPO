println("------------------------------------")
println("|              FockMPS             |")
println("------------------------------------")

@testset "FockMPS: arithmetic and canonicalize" begin
	L = 6
	D = 6
	tol = 1.0e-7
	for T in (Float64, ComplexF64)
		psi = randomfockmps(T, L, D=D)
		@test scalartype(psi) == T
		@test space_l(psi) == 1
		@test space_r(psi) == 1

		@test bond_dimension(psi) <= D
		psi1 = leftorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = tol
		@test distance(psi, psi1) / norm(psi) < tol

		psi1 = rightorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = tol
		@test distance(psi, psi1) / norm(psi) < tol

		psi1 = leftorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=true))
		@test isleftcanonical(psi1)
		psi1 = rightorth!(deepcopy(psi), alg = Orthogonalize(SVD(), normalize=true))
		@test isrightcanonical(psi1)
		psi1 = canonicalize!(deepcopy(psi), alg = Orthogonalize(SVD(), normalize=true))
		@test iscanonical(psi1)
		@test norm(2 * psi1) ≈ 2
		@test norm(psi1 / 2) ≈ 0.5
		@test norm(psi1 - psi1) ≈ 0. atol = tol
		@test distance(psi, psi) ≈ 0. atol = tol

		psi1 = canonicalize!(deepcopy(psi), alg=Orthogonalize(trunc=NoTruncation(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = tol
		@test distance(psi, psi1) / norm(psi) < tol

	end
	
end
