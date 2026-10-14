#default settings
module Defaults
	# a nested module does not see the enclosing module's bindings, so the
	# names used below are imported explicitly
	using Z2Tensors: TruncationScheme, truncdimcutoff, truncrelerr
	using ExpExp: OverDeterminedProny
	import ..DMRG1, ..WII
	const maxiter = 100 # for DMRG iteration
	const D = 100 # default bond dimension
	const tolgauge = 1e-14 # for MPS truncation
	const tol = 1e-12 # for DMRG iteration
	const tollanczos = 1.0e-10 # for lanczos eigensolver
	const tolexp = 1.0e-8 # for local eigen in DMRG
	const verbosity = 1

	alg_trunc_if(; D::Int=D, ϵ::Real=tolgauge, add_back::Int=0) = truncdimcutoff(D=D, ϵ=ϵ, add_back=add_back)
	alg_trunc_sys(; ϵ::Real=tolgauge) = truncrelerr(ϵ)
	alg_evo() = WII()
	alg_mult(; trunc::TruncationScheme=alg_trunc_if()) = DMRG1(trunc=trunc)
	alg_expan(; n::Int=15, tol::Real=1.0e-4, verbosity::Int=0) = OverDeterminedProny(n=n, tol=tol, verbosity=verbosity)
end
