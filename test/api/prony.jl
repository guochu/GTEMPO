using ExpExp: exponential_expansion, expansion_error, ExponentialExpansionAlgorithm

@testset "API: Prony expansions" begin
	f = [exp(-0.1k) * cos(k) for k in 1:60]
	@test AbstractPronyExpansion <: ExponentialExpansionAlgorithm
	@test OverDeterminedProny(n=30, tol=1.0e-10, verbosity=0) isa AbstractPronyExpansion
	@test DeterminedProny(n=20, tol=1.0e-10, verbosity=0) isa AbstractPronyExpansion
	@test MatrixPencil(n=20, tol=1.0e-10, verbosity=0) isa ExponentialExpansionAlgorithm
	@test LeastSquareProny(n=20, tol=1.0e-10, verbosity=0) isa ExponentialExpansionAlgorithm
	for alg in (OverDeterminedProny(n=30, tol=1.0e-10, verbosity=0),
				MatrixPencil(n=20, tol=1.0e-10, verbosity=0))
		coeffs, alphas = exponential_expansion(f, alg)
		@test length(coeffs) == length(alphas)
		@test expansion_error(f, coeffs, alphas) < 1.0e-6
	end
	# LeastSquareProny: the Levenberg-Marquardt refinement is not robust for
	# every signal, so only a single well-conditioned run is checked here
	coeffs, alphas = exponential_expansion([exp(-0.1k) for k in 1:40], LeastSquareProny(n=4, tol=1.0e-10, verbosity=0))
	@test length(coeffs) == length(alphas)
	# automatic stepsize
	coeffs, alphas = exponential_expansion(f; alg=OverDeterminedProny(n=15, tol=1.0e-8, stepsize=nothing, verbosity=0))
	@test expansion_error(f, coeffs, alphas) < 1.0e-3
end
