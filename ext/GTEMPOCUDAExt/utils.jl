
struct CuSVDCompression <: DMRGAlgorithm
	parent::SVDCompression
	CuSVDCompression(alg::SVDCompression) = new(alg)
end
Cu(alg::SVDCompression) = CuSVDCompression(alg)
CuSVDCompression(args...; kwargs...) = Cu(SVDCompression(args...; kwargs...))

abstract type CuDMRGMultAlgorithm <: DMRGAlgorithm end

struct CuDMRGMult1 <: CuDMRGMultAlgorithm
    parent::DMRGMult1
	CuDMRGMult1(alg::DMRGMult1) = new(alg)
end
Cu(alg::DMRGMult1) = CuDMRGMult1(alg)
CuDMRGMult1(args...; kwargs...) = Cu(DMRGMult1(args...; kwargs...))


CuAlgs = Union{CuSVDCompression, CuDMRGMult1}
function Base.getproperty(x::CuAlgs, s::Symbol)
	return s === :parent ? getfield(x, s) : getproperty(getfield(x, :parent), s)
end





function Z2TensorsCUDAExt.tocu(t::GrassmannTensorMap)
	return GrassmannTensorMap(tocu(t.data))
end
function Z2TensorsCUDAExt.fromcu(t::GrassmannTensorMap)
	return GrassmannTensorMap(fromcu(t.data))
end

function Z2TensorsCUDAExt.tocu(m::GrassmannTransferMatrix)
    states = map(x->map(i->tocu(i).data, x), m.states)
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
		u, s, v, err = stable_tsvd(GrassmannTensorMap(psii), (1,), (2, 3), trunc=trunc)
		psi[i] = get_data(fromcu(permute(v, (1,2), (3,))))
		nr = _renormalize!(psi, get_data(s), normalize)
		rerror = sqrt(err * err / (nr * nr + err * err))
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", rerror)
		u2 = u * s
		# nl = norm(u2)
		# (nl ≈ zero(nl)) && @warn "norm of GrassmannMPS is zero"
		# _rescaling!(psi, nl)
		# u2 = rmul!(u2, 1/nl)
		@tensor tmp[-1 -2; -3] := GrassmannTensorMap(tocu(psi[i-1]))[-1, -2, 1] * u2[1, -3]
		psii = get_data(tmp)
		psi.s[i] = fromcu(get_data(s))
		maxerr = max(maxerr, err)
	end
	psi[1] = fromcu(psii)
	(verbosity > 0) && println("Max SVD truncerror in rightorth: ", maxerr)
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

