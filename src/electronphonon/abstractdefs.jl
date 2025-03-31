abstract type Dense1DTN{T<:Number} end

const ValidIndices = Union{Integer,AbstractRange{Int64}, Colon}

TK.scalartype(::Type{<:Dense1DTN{T}}) where {T} = T
TK.scalartype(t::Dense1DTN) = scalartype(typeof(t))
Base.length(t::Dense1DTN) = length(t.data)
Base.getindex(t::Dense1DTN, i::ValidIndices) = getindex(t.data, i)
Base.setindex!(t::Dense1DTN, v, i::ValidIndices) = setindex!(t.data, v, i)
Base.firstindex(x::Dense1DTN) = firstindex(x.data)
Base.lastindex(t::Dense1DTN) = lastindex(t.data)
Base.isempty(t::Dense1DTN) = isempty(t.data)


const DenseMPSTensor{T} = AbstractArray{T, 3} where {T<:Number}
const DenseMPOTensor{T} = AbstractArray{T, 4} where {T<:Number}

DMRG.space_l(m::DenseMPSTensor) = size(m, 1)
DMRG.space_r(m::DenseMPSTensor) = size(m, 3)
DMRG.space_l(m::DenseMPOTensor) = size(m, 1)
DMRG.space_r(m::DenseMPOTensor) = size(m, 3)
DMRG.space_l(t::Dense1DTN) = size(t[1], 1)
DMRG.space_r(t::Dense1DTN) = size(t[end], 3)

DMRG.r_RR(psiA::Dense1DTN, psiB::Dense1DTN) = _eye(promote_type(scalartype(psiA), scalartype(psiB)), space_r(psiA), space_r(psiB))
DMRG.r_RR(psi::Dense1DTN) = r_RR(psi, psi)
DMRG.l_LL(psiA::Dense1DTN, psiB::Dense1DTN) = _eye(promote_type(scalartype(psiA), scalartype(psiB)), space_l(psiA), space_l(psiB))
DMRG.l_LL(psi::Dense1DTN) = l_LL(psi, psi)


DMRG.bond_dimension(psi::Dense1DTN, bond::Int) = begin
	((bond >= 1) && (bond <= length(psi))) || throw(BoundsError())
	space_r(psi[bond])
end 
DMRG.bond_dimensions(psi::Dense1DTN) = [bond_dimension(psi, i) for i in 1:length(psi)]
DMRG.bond_dimension(psi::Dense1DTN) = maximum(bond_dimensions(psi))


function DMRG.isleftcanonical(psij::DenseMPSTensor; kwargs...)
	@tensor r[-1, -2] := conj(psij[1,2,-1]) * psij[1,2,-2]
	return isapprox(r, one(r); kwargs...) 
end
function DMRG.isrightcanonical(psij::DenseMPSTensor; kwargs...)
	@tensor r[-1, -2] := conj(psij[-1,1,2]) * psij[-2,1,2]
	return isapprox(r, one(r); kwargs...) 
end
function DMRG.isleftcanonical(psij::DenseMPOTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[1,2,-1,3]) * psij[1,2,-2,3]
	return isapprox(r, one(r); kwargs...) 
end
function DMRG.isrightcanonical(psij::DenseMPOTensor; kwargs...)
	@tensor r[-1; -2] := conj(psij[-1,1,2,3]) * psij[-2,1,2,3]
	return isapprox(r, one(r); kwargs...) 
end