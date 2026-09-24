
"""
	FockMPO{T<:Number}

Fock-space finite matrix product operator: a lightweight wrapper of the
FiniteMPSAlgorithms `MPO` (field `.parent`), a chain of rank-4 site tensors.

Site tensor convention (payload `MPO`):
    o 
    |
    2
o-1   3-i
	4
	|
	i
i.e. dimension 1: left auxiliary (bond) index, dimension 2: output physical index,
dimension 3: right auxiliary (bond) index, dimension 4: input physical index.
The left and right boundaries are always vacuum (dimension 1).
`.data` delegates to the payload's site-tensor vector.
"""
struct FockMPO{T<:Number} <: Dense1DTN{T}
	parent::FMA.MPO{T}
end

function Base.getproperty(h::FockMPO, s::Symbol)
	s === :parent && return getfield(h, :parent)
	s === :data && return getfield(h, :parent).data
	throw(ArgumentError("FockMPO has no property $s"))
end
Base.propertynames(::FockMPO) = (:parent, :data)

FockMPO(data::AbstractVector{<:DenseMPOTensor{T}}) where {T<:Number} = FockMPO(FMA.MPO(data))

Base.copy(h::FockMPO) = FockMPO(copy(h.parent))

# apply the operator to a FockMPS (exact, no truncation)
Base.:*(h::FockMPO, psi::FockMPS) = FockMPS(h.parent * psi.parent)
