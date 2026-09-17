abstract type MPSAlgorithm end
abstract type DMRGAlgorithm <: MPSAlgorithm end

# truncation schemes carrying an explicit bond dimension D, required by the
# initial guesses (`:svd`/`:rand`/`:pre`) of the DMRG1/DMRG2 iterative multiplication
const TruncationWithD = Union{TruncationDimension, TruncateDimCutoff}

struct SVDCompression{T<:TruncationScheme} <: DMRGAlgorithm
	trunc::T
	verbosity::Int
end
SVDCompression(trunc::TruncationScheme; verbosity::Int=0) = SVDCompression(trunc, verbosity)
SVDCompression(; trunc::TruncationScheme=truncdimcutoff(D=Defaults.D, ϵ=Defaults.tol, add_back=0), verbosity::Int=0) = SVDCompression(trunc, verbosity)
Base.similar(x::SVDCompression; trunc::TruncationScheme=x.trunc, verbosity::Int=x.verbosity) = SVDCompression(trunc, verbosity=verbosity)

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
	trunc::TruncationWithD
	maxiter::Int
	tol::Float64
	initguess::Symbol
	verbosity::Int
	callback::Function
end
function DMRG1(trunc::TruncationWithD; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0, callback::Function=Returns(nothing))
	(initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
	return DMRG1(trunc, maxiter, tol, initguess, verbosity, callback)
end
DMRG1(; trunc::TruncationWithD=DefaultITruncation, kwargs...) = DMRG1(trunc; kwargs...)
Base.similar(x::DMRG1; trunc::TruncationWithD=x.trunc, maxiter::Int=x.maxiter, tol::Float64=x.tol, initguess::Symbol=x.initguess, verbosity::Int=x.verbosity, callback=x.callback) = DMRG1(
			trunc=trunc, maxiter=maxiter, tol=tol, initguess=initguess, verbosity=verbosity, callback=callback)

struct DMRG2 <: DMRGAlgorithm
	trunc::TruncationWithD
	maxiter::Int
	tol::Float64
	initguess::Symbol
	verbosity::Int
end

function DMRG2(trunc::TruncationWithD; maxiter::Int=5, tol::Float64=1.0e-12, initguess::Symbol=:svd, verbosity::Int=0)
	(initguess in AllowedInitGuesses) || throw(ArgumentError("initguess must be one of $(AllowedInitGuesses)"))
	return DMRG2(trunc, maxiter, tol, initguess, verbosity)
end
DMRG2(; trunc::TruncationWithD=DefaultITruncation, kwargs...) = DMRG2(trunc; kwargs...)
Base.similar(x::DMRG2; trunc::TruncationWithD=x.trunc, maxiter::Int=x.maxiter, tol::Float64=x.tol, initguess::Symbol=x.initguess, verbosity::Int=x.verbosity) = DMRG2(
			trunc=trunc, maxiter=maxiter, tol=tol, initguess=initguess, verbosity=verbosity)
