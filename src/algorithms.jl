abstract type MPSAlgorithm end
abstract type DMRGAlgorithm <: MPSAlgorithm end



struct SVDCompression <: DMRGAlgorithm
	D::Int
	tol::Float64
	verbosity::Int
end
SVDCompression(; D::Int=Defaults.D, tol::Float64=Defaults.tol, verbosity::Int=0) = SVDCompression(D, tol, verbosity)

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

# orthogonalize mps to be left-canonical or right-canonical
abstract type MatrixProductOrthogonalAlgorithm  end

"""
	struct MatrixProductOrthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme}
"""
struct Orthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme} <: MatrixProductOrthogonalAlgorithm
	orth::A
	trunc::T
	normalize::Bool
	verbosity::Int
end
Orthogonalize(a::Union{QR, SVD}, trunc::TruncationScheme; normalize::Bool=false, verbosity::Int=0) = Orthogonalize(a, trunc, normalize, verbosity)
Orthogonalize(a::Union{QR, SVD}; trunc::TruncationScheme=TK.NoTruncation(), normalize::Bool=false, verbosity::Int=0) = Orthogonalize(a, trunc, normalize, verbosity)
Orthogonalize(; alg::Union{QR, SVD} = SVD(), trunc::TruncationScheme=TK.NoTruncation(), normalize::Bool=false, verbosity::Int=0) = Orthogonalize(alg, trunc, normalize, verbosity)

const AllowedInitGuesses = (:svd, :pre, :rand)

struct DMRG1 <: DMRGAlgorithm
	trunc::TruncationDimCutoff
	maxiter::Int
	tol::Float64
	initguess::Symbol
	verbosity::Int
	callback::Function
end
function DMRG1(trunc::TruncationDimCutoff; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0, callback::Function=Returns(nothing))
	(initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
	return DMRG1(trunc, maxiter, tol, initguess, verbosity, callback)
end
DMRG1(; trunc::TruncationDimCutoff=DefaultITruncation, kwargs...) = DMRG1(trunc; kwargs...)
Base.similar(x::DMRG1; trunc::TruncationDimCutoff=x.trunc, maxiter::Int=x.maxiter, tol::Float64=x.tol, initguess::Symbol=x.initguess, verbosity::Int=x.verbosity, callback=x.callback) = DMRG1(
			trunc=trunc, maxiter=maxiter, tol=tol, initguess=initguess, verbosity=verbosity, callback=callback)

struct DMRG2 <: DMRGAlgorithm
	trunc::TruncationDimCutoff
	maxiter::Int
	tol::Float64
	initguess::Symbol
	verbosity::Int
end

function DMRG2(trunc::TruncationDimCutoff; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0)
	(initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
	return DMRG2(trunc, maxiter, tol, initguess, verbosity)
end
DMRG2(; trunc::TruncationDimCutoff=DefaultITruncation, kwargs...) = DMRG2(trunc; kwargs...)
Base.similar(x::DMRG2; trunc::TruncationDimCutoff=x.trunc, maxiter::Int=x.maxiter, tol::Float64=x.tol, initguess::Symbol=x.initguess, verbosity::Int=x.verbosity) = DMRG2(
			trunc=trunc, maxiter=maxiter, tol=tol, initguess=initguess, verbosity=verbosity)

function Base.getproperty(x::DMRGAlgorithm, s::Symbol)
	if s == :D
		return x.trunc.D
	elseif s == :ϵ
		return x.trunc.ϵ
	else
		getfield(x, s)
	end
end
