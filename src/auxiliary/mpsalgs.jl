abstract type MPSAlgorithm end
abstract type DMRGAlgorithm <: MPSAlgorithm end



@with_kw struct SVDCompression <: DMRGAlgorithm
	D::Int = Defaults.D
	tol::Float64 = Defaults.tol
	verbosity::Int = 0
end

# @with_kw struct Deparallelise <: DMRGAlgorithm
# 	tol::Float64 = DeparalleliseTol
# 	verbosity::Int = Defaults.verbosity
# end

SVDCompression(trunc::TruncationDimCutoff; verbosity::Int=0) = SVDCompression(D=trunc.D, tol=trunc.ϵ, verbosity=verbosity)
Base.similar(x::SVDCompression; D::Int=x.D, tol::Float64=x.tol, verbosity::Int=x.verbosity) = SVDCompression(D=D, tol=tol, verbosity=verbosity)

function Base.getproperty(x::SVDCompression, s::Symbol)
	if s == :trunc
		return get_trunc(x)
	elseif s == :ϵ
		return x.tol
	else
		getfield(x, s)
	end
end

get_trunc(alg::SVDCompression) = truncdimcutoff(D=alg.D, ϵ=alg.tol, add_back=0)

# compress!(h::MPO, alg::SVDCompression) = canonicalize!(h, alg=Orthogonalize(SVD(), get_trunc(alg); normalize=false))
# compress!(h::MPO, alg::Deparallelise) = deparallel!(h, tol=alg.tol, verbosity=alg.verbosity)
# compress!(h::MPO; alg::DMRGAlgorithm = Deparallelise()) = compress!(h, alg)
# compress!(psi::MPS, alg::SVDCompression) = canonicalize!(psi, alg=Orthogonalize(trunc=get_trunc(alg), normalize=false))
