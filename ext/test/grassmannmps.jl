println("------------------------------------")
println("|            GrassmannMPS          |")
println("------------------------------------")

@testset "GrassmannMPS: multiplications" begin
	L = 6
	chi = 20
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	alg1 = CuSVDCompression(trunc)
	alg2 = CuDMRGMult1(trunc, initguess=:svd)
	alg3 = CuDMRGMult1(trunc, initguess=:rand, maxiter=10)
	alg4 = CuDMRGMult1(trunc, initguess=:pre, maxiter=10)
	algs = [alg1, alg2, alg3, alg4]
	tol = 1.0e-7
	for T in (Float64, ComplexF64)
		psi1 = randomgmps(T, L, D=4)
		psi2 = randomgmps(T, L, D=4)

		psi3 = psi1 * psi2
		_n = norm(psi3)
		for alg in algs
			psi4 = mult(psi1, psi2, alg)
			@test distance(psi3, psi4) / _n < tol
		end

		canonicalize!(psi1)
		canonicalize!(psi2)
		psi5 = psi1 * psi2
		@test distance(psi3, psi5) / _n < tol

		for alg in algs
			psi4 = mult(psi1, psi2, alg)
			@test distance(psi3, psi4) / _n < tol
		end
	end
end

