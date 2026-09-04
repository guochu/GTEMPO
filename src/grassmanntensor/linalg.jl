#---------------------------------------------------------------
# Grassmann (fermionic-sign) linear algebra helpers on plain parity
# TensorMaps: the fermionic signs enter through `f_permute`
# (see grassmanntensor.jl); the decompositions themselves are the
# ordinary (bosonic) TensorMap ones.
#---------------------------------------------------------------

# fermionic permute
g_permute(t::AbstractParityTensorMap, p::Index2Tuple; copy::Bool=false) = f_permute(t, p; copy=copy)
g_permute(t::AbstractParityTensorMap, p1::IndexTuple, p2::IndexTuple; copy::Bool=false) = f_permute(t, (p1, p2); copy=copy)

# stable tsvd after a fermionic reordering
function g_stable_tsvd(t::AbstractParityTensorMap, p1::IndexTuple, p2::IndexTuple; kwargs...)
	return stable_tsvd!(g_permute(t, p1, p2; copy=true); kwargs...)
end

# rightorth after a fermionic reordering
g_rightorth(t::AbstractParityTensorMap, p1::IndexTuple, p2::IndexTuple; kwargs...) = TK.rightorth!(g_permute(t, p1, p2; copy=true); kwargs...)
