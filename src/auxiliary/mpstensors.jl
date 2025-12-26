const MPSTensor{S} = AbstractTensorMap{<:Number, S, 2, 1} where {S<:ElementarySpace}
const MPOTensor{S} = AbstractTensorMap{<:Number, S, 2, 2} where {S<:ElementarySpace}
# The bond tensor are just the singlur vector but has to be stored as a general matrix 
# since TensorKit does not specialize for Diagonal Matrices
const MPSBondTensor{S} = AbstractTensorMap{<:Number, S, 1, 1} where {S<:ElementarySpace}
const SiteOperator{S} = Union{MPOTensor{S}, MPSBondTensor{S}}


mpstensortype(::Type{S}, ::Type{T}) where {S <: ElementarySpace, T} = tensormaptype(S, 2, 1, T)
mpotensortype(::Type{S}, ::Type{T}) where {S <: ElementarySpace, T} = tensormaptype(S, 2, 2, T)
function bondtensortype(::Type{S}, ::Type{TorA}) where {S <: ElementarySpace, TorA<:Union{Number, DenseVector}} 
    if TorA <: Number
        return DiagonalTensorMap{TorA,S,Vector{TorA}}
    elseif TorA <: DenseVector
        return DiagonalTensorMap{scalartype(TorA),S,TorA}
    else
        throw(ArgumentError("argument $TorA should specify a scalar type (`<:Number`) or a storage type `<:DenseVector{<:Number}`"))
    end
end


space_l(a::MPSBondTensor) = space(a, 1)
space_r(a::MPSBondTensor) = space(a, 2)
space_l(a::MPSTensor) = space(a, 1)
space_r(a::MPSTensor) = space(a, 3)
space_l(a::MPOTensor) = space(a, 1)
space_r(a::MPOTensor) = space(a, 3)
physical_space(a::MPSTensor) = space(a, 2)
ophysical_space(a::MPOTensor) = space(a, 2)
iphysical_space(a::MPOTensor) = space(a, 4)
ophysical_space(a::MPSBondTensor) = space(a, 1)
iphysical_space(a::MPSBondTensor) = space(a, 2)
physical_space(a::SiteOperator) = ophysical_space(a)

"""
	r_RR(a::MPS, b::MPS)

Notice the convention!!!
a is bra, b is ket, ^
					a
					-
					b
					v
for r_RR b is codomain, namely 
		----2
r_RR =
		----1
for l_LL a is codomain, namely
		----1
l_LL =
		----2
"""
# r_RR(a::MPSTensor, b::MPSTensor) = loose_isometry(Matrix{promote_type(scalartype(a), scalartype(b))}, space_r(b)', space_r(a)')
# l_LL(a::MPSTensor, b::MPSTensor) = loose_isometry(Matrix{promote_type(scalartype(a), scalartype(b))}, space_l(a), space_l(b))
function l_LL(a::MPSTensor, b::MPSTensor)
	@assert isoneunit(space_l(a)) && isoneunit(space_l(b))
	ones(promote_type(scalartype(a), scalartype(b)), space_l(a), space_l(b))
end

isoneunit(s::ElementarySpace) = isdual(s) ? dual(s) == oneunit(s) : s == oneunit(s) 


# check if mps tensor is canonical
function isleftcanonical(psij::MPSTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[1,2,-1]) * psij[1,2,-2]
	return isapprox(r, one(r); kwargs...) 
end
function isleftcanonical_r(psij::MPSTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[1,2,-1]) * psij[1,2,-2]
	for (c, b) in blocks(r)
		_one = lmul!( 1/(dim(c) * size(b, 1)), one(b))
		isapprox(b, _one; kwargs...) || return false
	end
	return true
end
function isrightcanonical(psij::MPSTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[-1,1,2]) * psij[-2,1,2]
	return isapprox(r, one(r); kwargs...) 
end

# check if mpo tensor is canonical
function isleftcanonical(psij::MPOTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[1,2,-1,3]) * psij[1,2,-2,3]
	return isapprox(r, one(r); kwargs...) 
end
function isleftcanonical_r(psij::MPOTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[1,2,-1,3]) * psij[1,2,-2,3]
	for (c, b) in blocks(r)
		_one = lmul!( 1/(dim(c) * size(b, 1)), one(b))
		isapprox(b, _one; kwargs...) || return false
	end
	return true
end
function isrightcanonical(psij::MPOTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[-1,1,2,3]) * psij[-2,1,2,3]
	return isapprox(r, one(r); kwargs...) 
end

isstrict(s::ElementarySpace) = isoneunit(s) || (FusionStyle(sectortype(s)) isa UniqueFusion)




"""
	updateright(hold::MPSBondTensor, mpsAj::MPSTensor{S}, mpsBj::MPSTensor{S}) where {S<:ElementarySpace}
	update storage from right to left for overlap of mps
"""
function updateright(hold::MPSBondTensor, mpsAj::MPSTensor, mpsBj::MPSTensor) 
	@tensor hnew[1;5] := mpsBj[1,2,3] * hold[3,4] * conj(mpsAj[5,2,4])
end

"""
	updateleft(hold::MPSBondTensor, mpsAj::MPSTensor{S}, mpsBj::MPSTensor{S}) where {S<:ElementarySpace}
	update storage from left to right for overlap of mps
"""
function updateleft(hold::MPSBondTensor, mpsAj::MPSTensor, mpsBj::MPSTensor) 
	@tensor m2[-1 -2; -3] := conj(mpsAj[1, -2, -3]) * hold[1, -1]
	@tensor hnew[-1; -2] := m2[1,2,-1] * mpsBj[1,2,-2]
end

