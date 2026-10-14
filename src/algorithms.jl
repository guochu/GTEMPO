abstract type MPSAlgorithm end
abstract type DMRGAlgorithm <: MPSAlgorithm end

# truncation schemes carrying an explicit bond dimension D, required by the
# initial guesses (`:svd`/`:rand`/`:pre`) of the DMRG1/DMRG2 iterative multiplication
const TruncationWithD = Union{TruncationDimension, TruncateDimCutoff}

@kwdef struct SVDCompression{T<:TruncationScheme} <: DMRGAlgorithm
	trunc::T = Defaults.alg_trunc_if(ϵ=Defaults.tol)
	verbosity::Int = 0
end

# compress!(h::MPO, alg::SVDCompression) = canonicalize!(h, alg=Orthogonalize(SVD(), get_trunc(alg); normalize=false))
# compress!(h::MPO, alg::Deparallelise) = deparallel!(h, tol=alg.tol, verbosity=alg.verbosity)
# compress!(h::MPO; alg::DMRGAlgorithm = Deparallelise()) = compress!(h, alg)
# compress!(psi::MPS, alg::SVDCompression) = canonicalize!(psi, alg=Orthogonalize(trunc=get_trunc(alg), normalize=false))

# orthogonalize mps to be left-canonical or right-canonical
abstract type MatrixProductOrthogonalAlgorithm  end

"""
	struct MatrixProductOrthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme}
"""
@kwdef struct Orthogonalize{A<:Union{QR, SVD}, T<:TruncationScheme} <: MatrixProductOrthogonalAlgorithm
	orth::A = SVD()
	trunc::T = Defaults.alg_trunc_sys()
	normalize::Bool = false
	verbosity::Int = 0
end

@kwdef struct DMRG1 <: DMRGAlgorithm
	trunc::TruncationWithD = Defaults.alg_trunc_if()
	maxiter::Int = 5
	tol::Float64 = 1.0e-12
	initguess::Symbol = :svd
	verbosity::Int = 0
end

@kwdef struct DMRG2 <: DMRGAlgorithm
	trunc::TruncationWithD = Defaults.alg_trunc_if()
	maxiter::Int = 5
	tol::Float64 = 1.0e-12
	initguess::Symbol = :svd
	verbosity::Int = 0
end
