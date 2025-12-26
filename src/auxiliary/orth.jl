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
