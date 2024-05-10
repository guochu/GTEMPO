Base.:*(x::GrassmannTensorMap, y::GrassmannTensorMap) = GrassmannTensorMap(x.data * y.data)

function DMRG.stable_tsvd!(t::GrassmannTensorMap; kwargs...)
	u, s, v, err = stable_tsvd!(t.data; kwargs...)
	return GrassmannTensorMap(u), GrassmannTensorMap(s), GrassmannTensorMap(v), err
end 
DMRG.stable_tsvd(t::GrassmannTensorMap, p1::IndexTuple, p2::IndexTuple; kwargs...) = stable_tsvd!(permute(t, p1, p2, copy=true); kwargs...)