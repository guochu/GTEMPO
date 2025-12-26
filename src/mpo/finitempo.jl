# another MPO object is needed to support non-number conserving operators, such a

"""
	MPO{A <: MPOTensor}
Finite Matrix Product Operator which stores a chain of rank-4 site tensors.
"""
struct MPO{A <: MPOTensor} <: AbstractFiniteMPO{A}
	data::Vector{A}

"""
	MPO{A}(mpotensors::Vector)
Constructor entrance for MPO, which only supports strictly quantum number conserving operators

Site tensor convention:
i mean in arrow, o means out arrow
    o 
    |
    2
o-1   3-i
	4
	|
	i
The left and right boundaries are always vacuum.
The case that the right boundary is not vacuum corresponds to operators which do 
not conserve quantum number, such as a†.
"""
function MPO(data::Vector{A}) where {A<:MPOTensor}
	@assert !isempty(data)
	check_mpo_spaces(data)
	return new{A}(data)
end
end

function check_mpo_spaces(mpotensors::AbstractVector)
	# all(check_mpotensor_dir, mpotensors) || throw(SpaceMismatch())
	for i in 1:length(mpotensors)-1
		(space_r(mpotensors[i]) == space_l(mpotensors[i+1])') || throw(SpaceMismatch())
	end
	isoneunit(space_l(mpotensors[1])) || throw(SpaceMismatch("space_l of the left boundary tensor should be vacuum by convention"))
end

storage(a::MPO) = a.data
Base.length(a::MPO) = length(storage(a))
Base.isempty(a::MPO) = isempty(storage(a))
Base.getindex(a::MPO, i::Int) = getindex(storage(a), i)
Base.firstindex(a::MPO) = firstindex(storage(a))
Base.lastindex(a::MPO) = lastindex(storage(a))


function Base.setindex!(h::MPO, v::MPOTensor, i::Int)
	# check_mpotensor_dir(v) || throw(SpaceMismatch())
	if i == 1
		isoneunit(space_l(v)) || throw(SpaceMismatch("space_l of the left boundary tensor should be vacuum by convention"))
	end
	return setindex!(h.data, v, i)
end 
Base.copy(h::MPO) = MPO(copy(h.data))

