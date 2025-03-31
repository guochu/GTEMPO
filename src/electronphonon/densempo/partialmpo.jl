# another MPO object is needed to support non-number conserving operators, such a

"""
	PartialDenseMPO{A <: MPOTensor}
Finite Matrix Product Operator which stores a chain of rank-4 site tensors.
"""
struct PartialDenseMPO{T<:Number} 
	data::Vector{Array{T, 4}}
	positions::Vector{Int}

"""
	PartialDenseMPO{A}(mpotensors::Vector, positions::Vector{Int})

This is mainly used for computing observables
"""
function PartialDenseMPO(data::AbstractVector{<:DenseMPOTensor{T}}, positions::AbstractVector{Int}) where {T<:Number}
	@assert !isempty(data)
	@assert length(data) == length(positions)
	# @assert positions[1] > 0 # allow periodic condition
	_check_mpo_space(data)
	for i in 1:length(positions)-1
		(positions[i] < positions[i+1]) || throw(ArgumentError("positions should be ordered from small to large"))
	end
	return new{T}(convert(Vector{Array{T, 4}}, data), convert(Vector{Int}, positions))
end
end

PartialDenseMPO(positions::AbstractVector{Int}, data::AbstractVector{<:DenseMPOTensor}) = PartialDenseMPO(data, positions)

TK.scalartype(::Type{PartialDenseMPO{T}}) where {T<:Number} = T

storage(a::PartialDenseMPO) = a.data
DMRG.positions(a::PartialDenseMPO) = a.positions
Base.length(a::PartialDenseMPO) = length(storage(a))
Base.isempty(a::PartialDenseMPO) = isempty(storage(a))
Base.getindex(a::PartialDenseMPO, i::Int) = getindex(storage(a), i)
Base.firstindex(a::PartialDenseMPO) = firstindex(storage(a))
Base.lastindex(a::PartialDenseMPO) = lastindex(storage(a))

DMRG.space_l(state::PartialDenseMPO) = space_l(state[1])
DMRG.space_r(state::PartialDenseMPO) = space_r(state[end])

function Base.setindex!(h::PartialDenseMPO, v::DenseMPOTensor, i::Int)
	# check_mpotensor_dir(v) || throw(SpaceMismatch())
	if i == 1
		(space_l(v) == 1) || throw(SpaceMismatch("space_l of the left boundary tensor should be vacuum by convention"))
	end
	return setindex!(h.data, v, i)
end 
Base.copy(h::PartialDenseMPO) = PartialDenseMPO(copy(storage(h)), copy(positions(h)))

function Base.complex(psi::PartialDenseMPO)
	if scalartype(psi) <: Real
		data = [complex(item) for item in psi.data]
		return PartialDenseMPO(data, positions(psi))
	end
	return psi
end

TK.id(m::PartialDenseMPO) = PartialDenseMPO([reshape(_eye(scalartype(m), 2, 2), 1,2,1,2) for item in m.data], positions(m))
