@testset "API: mult" begin
	# two random states
	a = randomgmps(6, D=4)
	b = randomgmps(6, D=4)
	exact = a * b
	relerr(z) = distance(z, exact) / norm(exact)

	# with a truncation large enough, mult(a, b) reproduces the exact product
	truncbig = truncdimcutoff(D=256, ϵ=1.0e-12)
	z = mult(a, b, trunc=truncbig)[1]
	@test relerr(z) < 1.0e-6
	@test iscanonical(z)

	z, maxerr = mult(a, b, SVDCompression(truncbig))
	@test relerr(z) < 1.0e-6
	@test iscanonical(z)
	@test maxerr < 1.0e-6

	z, info = mult(a, b, DMRG1(truncbig, maxiter=30))
	@test relerr(z) < 1.0e-6
	@test iscanonical(z)
	@test info.converged

	z, info = mult(a, b, DMRG2(truncbig, maxiter=30))
	@test relerr(z) < 1.0e-6
	@test iscanonical(z)
	@test info.converged

	# mult! stores the result in the first argument
	x = copy(a)
	x2, info = mult!(x, b, DMRG1(truncbig, maxiter=30))
	@test x2 === x
	@test relerr(x) < 1.0e-6

	# with an active truncation the mixed-canonical form is not exact, only
	# check that the tensors stay right-canonical
	for alg in (SVDCompression(truncdimcutoff(D=2, ϵ=1.0e-8)),
				DMRG1(truncdimcutoff(D=2, ϵ=1.0e-8), verbosity=-1),
				DMRG2(truncdimcutoff(D=2, ϵ=1.0e-8), verbosity=-1))
		z = mult(a, b, alg)[1]
		@test isrightcanonical(z)
	end
end
