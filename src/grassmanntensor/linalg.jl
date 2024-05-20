Base.:*(x::GrassmannTensorMap, y::GrassmannTensorMap) = GrassmannTensorMap(x.data * y.data)
TK.lmul!(x::Number, t::GrassmannTensorMap) = GrassmannTensorMap(lmul!(x, t.data))
TK.rmul!(t::GrassmannTensorMap, x::Number) = GrassmannTensorMap(rmul!(t.data, x))


function DMRG.stable_tsvd!(t::GrassmannTensorMap; kwargs...)
	u, s, v, err = stable_tsvd!(t.data; kwargs...)
	return GrassmannTensorMap(u), GrassmannTensorMap(s), GrassmannTensorMap(v), err
end 
DMRG.stable_tsvd(t::GrassmannTensorMap, p1::IndexTuple, p2::IndexTuple; kwargs...) = stable_tsvd!(permute(t, p1, p2, copy=true); kwargs...)
function TK.leftorth!(t::GrassmannTensorMap; kwargs...)
	q, r = leftorth!(t.data; kwargs...)
	return GrassmannTensorMap(q), GrassmannTensorMap(r)
end
function TK.rightorth!(t::GrassmannTensorMap; kwargs...)
	l, q = rightorth!(t.data; kwargs...)
	return GrassmannTensorMap(l), GrassmannTensorMap(q)
end
TK.rightorth(t::GrassmannTensorMap, p1::IndexTuple, p2::IndexTuple; kwargs...) = rightorth!(permute(t, p1, p2, copy=true); kwargs...)