# default settings
module Defaults
	# a nested module does not see the enclosing module's bindings, so the
	# names used below are imported explicitly
	using GrassmannMPS: TruncationStrategy, truncdimcutoff, truncrelerr
	using ExpExp: OverDeterminedProny
	import ..SVDCompression

	const D = 128       # default bond dimension
	const rtol = 1e-12  # default relative truncation tolerance
	const tolgauge = 1e-14
	verbosity = 1

	# 算法缺省（backup/GTEMPO defaults.jl 移植；alg_mult 适配 GM 后端
	# SVDCompression，TruncationScheme → TruncationStrategy 随 MAK 体系更名）：
	alg_trunc_if(; D::Int=D, ϵ::Real=tolgauge, add_back::Int=0) = truncdimcutoff(D=D, ϵ=ϵ, add_back=add_back)
	alg_trunc_sys(; ϵ::Real=tolgauge) = truncrelerr(ϵ)
	alg_mult(; trunc::TruncationStrategy=alg_trunc_if()) = SVDCompression(trunc=trunc)
	alg_expan(; n::Int=15, tol::Real=1.0e-4, verbosity::Int=0) = OverDeterminedProny(n=n, tol=tol, verbosity=verbosity)
end
