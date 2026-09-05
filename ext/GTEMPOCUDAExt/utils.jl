
struct CuSVDCompression <: DMRGAlgorithm
	parent::SVDCompression
	CuSVDCompression(alg::SVDCompression) = new(alg)
end
Cu(alg::SVDCompression) = CuSVDCompression(alg)
CuSVDCompression(args...; kwargs...) = Cu(SVDCompression(args...; kwargs...))

abstract type CuDMRGMultAlgorithm <: DMRGAlgorithm end

struct CuDMRG1 <: CuDMRGMultAlgorithm
    parent::DMRG1
	CuDMRG1(alg::DMRG1) = new(alg)
end
Cu(alg::DMRG1) = CuDMRG1(alg)
CuDMRG1(args...; kwargs...) = Cu(DMRG1(args...; kwargs...))


CuAlgs = Union{CuSVDCompression, CuDMRG1}
function Base.getproperty(x::CuAlgs, s::Symbol)
	return s === :parent ? getfield(x, s) : getproperty(getfield(x, :parent), s)
end





function Z2TensorsCUDAExt.tocu(m::GrassmannTransferMatrix)
    states = map(x->tocu.(x), m.states)
    return GrassmannTransferMatrix(states, m.scaling)
end


const CuStridedView{T, N, A <: CuArray{T}} = StridedView{T, N, A}
function LinearAlgebra.axpy!(a::Number, X::CuStridedView{<:Number}, Y::CuStridedView{T,N}) where {T<:Number,N}
	TO.tensoradd!(Y, X, (ntuple(i->i,N,),()), false, a, 1)
end
function LinearAlgebra.lmul!(a::Number, X::CuStridedView{T,N}) where {T<:Number,N}
	TO.tensoradd!(X, X, (ntuple(i->i,N,),()), false, a, 0)
end





function _cu_rightorth!(psi::GrassmannMPS, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	L = length(psi)
	maxerr = 0.
	psii = tocu(psi[end])
	for i in L:-1:2
		# single-site operation: fermionic twist + restore cancel exactly,
		# so the plain bosonic tsvd on the permuted tensor is identical
		u, s, v, err = tsvd(permute(psii, (1,), (2, 3); copy=true); alg=SDD(), trunc=trunc)
		psi[i] = fromcu(permute(v, (1, 2), (3,); copy=true))
		nr = _renormalize!(psi, s, normalize)
		rerror = sqrt(err * err / (nr * nr + err * err))
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", rerror)
		u2 = u * s
		# nl = norm(u2)
		# (nl ≈ zero(nl)) && @warn "norm of GrassmannMPS is zero"
		# _rescaling!(psi, nl)
		# u2 = rmul!(u2, 1/nl)
		@grassmann tmp[-1 -2; -3] := tocu(psi[i-1])[-1, -2, 1] * u2[1, -3]
		psii = tmp
		psi.s[i] = fromcu(s)
		maxerr = max(maxerr, err)
	end
	psi[1] = fromcu(psii)
	(verbosity > 0) && println("Max SVD truncerror in rightorth: ", maxerr)
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

