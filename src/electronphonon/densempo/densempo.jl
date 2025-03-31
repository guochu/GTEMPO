
"""
	DenseMPO{A <: MPOTensor}
Finite Matrix Product Operator which stores a chain of rank-4 site tensors.
"""
struct DenseMPO{T<:Number} <: Dense1DTN{T}
	data::Vector{Array{T, 4}}

"""
	DenseMPO{A}(mpotensors::Vector)
Constructor entrance for MPO, which only supports strictly quantum number conserving operators

site tensor convention:
i mean in arrow, o means out arrow
    o 
    |
    2
o-1   3-i
	4
	|
	i
The left and right boundaries are always vacuum.
The case that the right boundary is not vacuum corresponds to operators which do not conserve quantum number, 
such as a†, this case is implemented with another MPO object.
"""
function DenseMPO{T}(mpotensors::AbstractVector) where {T<:Number}
	_check_mpo_space(mpotensors)
	return new{T}(convert(Vector{Array{T, 4}}, mpotensors))
end

end
DenseMPO(data::AbstractVector{<:DenseMPOTensor{T}}) where {T <: Number} = MPO{T}(data)


# attributes
ophysical_dimensions(psi::DenseMPO) = [size(item, 2) for item in psi.data]
iphysical_dimensions(psi::DenseMPO) = [size(item, 4) for item in psi.data]


function _check_mpo_space(mpotensors::Vector)
	L = length(mpotensors)
	for i in 1:L-1
		(space_r(mpotensors[i]) == space_l(mpotensors[i+1])) || throw(DimensionMismatch())
	end
	# boundaries should be dimension 
	(space_l(mpotensors[1]) == 1) || throw(DimensionMismatch())
	(space_r(mpotensors[L]) == 1) || throw(DimensionMismatch())
	all(x->size(x, 2)==size(x,4)==2, mpotensors) || throw(ArgumentError("physical dimension of site tensors must be 2"))
	return true
end
