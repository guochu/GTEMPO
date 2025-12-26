# another MPO object is needed to support non-number conserving operators, such a

"""
	PartialMPO{A <: MPOTensor}
Finite Matrix Product Operator which stores a chain of rank-4 site tensors.
"""
struct PartialMPO{A <: MPOTensor} 
	data::Vector{A}
	positions::Vector{Int}

"""
	PartialMPO{A}(mpotensors::Vector, positions::Vector{Int})

This is mainly used for computing observables
"""
function PartialMPO(data::AbstractVector{A}, positions::AbstractVector{Int}) where {A<:MPOTensor}
	@assert !isempty(data)
	@assert length(data) == length(positions)
	# @assert positions[1] > 0 # allow periodic condition
	check_mpo_spaces(data)
	isoneunit(space_r(data[end])) || throw(ArgumentError("only strict PartialMPO allowed"))
	for i in 1:length(positions)-1
		(positions[i] < positions[i+1]) || throw(ArgumentError("positions should be ordered from small to large"))
	end
	return new{A}(convert(Vector{A}, data), convert(Vector{Int}, positions))
end
end

PartialMPO(positions::AbstractVector{Int}, data::AbstractVector{<:MPOTensor}) = PartialMPO(data, positions)

TK.scalartype(::Type{PartialMPO{A}}) where {A<:MPOTensor} = scalartype(A)
TK.spacetype(::Type{PartialMPO{A}}) where {A<:MPOTensor} = spacetype(A)
TK.spacetype(m::PartialMPO) = spacetype(typeof(m))


storage(a::PartialMPO) = a.data
positions(a::PartialMPO) = a.positions
Base.length(a::PartialMPO) = length(storage(a))
Base.isempty(a::PartialMPO) = isempty(storage(a))
Base.getindex(a::PartialMPO, i::Int) = getindex(storage(a), i)
Base.firstindex(a::PartialMPO) = firstindex(storage(a))
Base.lastindex(a::PartialMPO) = lastindex(storage(a))

space_l(state::PartialMPO) = space_l(state[1])
space_r(state::PartialMPO) = space_r(state[end])

function Base.setindex!(h::PartialMPO, v::MPOTensor, i::Int)
	# check_mpotensor_dir(v) || throw(SpaceMismatch())
	if i == 1
		isoneunit(space_l(v)) || throw(SpaceMismatch("space_l of the left boundary tensor should be vacuum by convention"))
	end
	return setindex!(h.data, v, i)
end 
Base.copy(h::PartialMPO) = PartialMPO(copy(storage(h)), copy(positions(h)))

function Base.complex(psi::PartialMPO)
	if scalartype(psi) <: Real
		data = [complex(item) for item in psi.data]
		return PartialMPO(data, positions(psi))
	end
	return psi
end


ophysical_space(a::PartialMPO, i::Int) = ophysical_space(a[i])
iphysical_space(a::PartialMPO, i::Int) = iphysical_space(a[i])
function physical_spaces(psi::PartialMPO)
	xs = ophysical_spaces(psi)
	(xs == adjoint.(iphysical_spaces(psi))) || throw(SpaceMismatch("i and o physical dimension mismatch."))
	return xs
end
left_virtualspace(a::PartialMPO, i::Int) = space_l(a[i])
right_virtualspace(a::PartialMPO, i::Int) = space_r(a[i])
left_virtualspaces(a::PartialMPO) = [left_virtualspace(a, i) for i in 1:length(a)]
right_virtualspaces(a::PartialMPO) = [right_virtualspace(a, i) for i in 1:length(a)]


ophysical_spaces(psi::PartialMPO) = [ophysical_space(psi[i]) for i in 1:length(psi)]
iphysical_spaces(psi::PartialMPO) = [iphysical_space(psi[i]) for i in 1:length(psi)]

