# store mpo tensor as a matrix of mpotensors instead of a single mpotensor
abstract type AbstractSparseMPOTensor{S} end
const SparseMPOTensorElement{M<:MPOTensor, T<:Number} = Union{M, T}

TK.spacetype(::Type{<:AbstractSparseMPOTensor{S}}) where S = S
TK.spacetype(x::AbstractSparseMPOTensor) = spacetype(typeof(x))


function sparsempotensoreltype(::Type{S}, ::Type{T}) where {S <: ElementarySpace, T}
	M = mpotensortype(S, T)
	return SparseMPOTensorElement{M, eltype(M)}
end

storage(m::AbstractSparseMPOTensor) = m.Os

Base.size(m::AbstractSparseMPOTensor) = size(storage(m))
Base.size(m::AbstractSparseMPOTensor, i::Int) = size(storage(m), i)

function Base.getindex(m::AbstractSparseMPOTensor, j::Int, k::Int) 
	T = scalartype(m)
	r = getindex(storage(m), j, k)
	if isa(r, T)
		imspace = m.leftspaces[j]
		omspace = m.rightspaces[k]
		pspace = m.pspace
		if r == zero(T)
			return zeros(T, imspace ⊗ pspace, omspace ⊗ pspace)
		else
			return r * isomorphism(T, imspace ⊗ pspace, omspace ⊗ pspace)
		end
	else
		return r
	end
end 
Base.lastindex(m::AbstractSparseMPOTensor) = lastindex(storage(m))
Base.lastindex(m::AbstractSparseMPOTensor, i::Int) = lastindex(storage(m), i)

Base.keys(x::AbstractSparseMPOTensor) = Iterators.filter(a->contains(x, a[1],a[2]),Iterators.product(1:size(x, 1),1:size(x, 2)))
opkeys(x::AbstractSparseMPOTensor) = Iterators.filter(a-> !isscal(x,a[1],a[2]),keys(x))
scalkeys(x::AbstractSparseMPOTensor) = Iterators.filter(a-> isscal(x,a[1],a[2]),keys(x))


"""
	isid(x::MPOTensor; kwargs...)
	isid(x::MPSBondTensor; kwargs...)

Check if given MPOTensor or MPSBondTensor is identity 
"""
function isid(x::MPOTensor; kwargs...)
    cod = space(x,1) ⊗ space(x,2);
    dom = space(x,3)' ⊗ space(x,4)';

    #would like to have an 'isisomorphic'
    for c in union(blocksectors(cod), blocksectors(dom))
        blockdim(cod, c) == blockdim(dom, c) || return false,0.0;
    end

    id = isomorphism(scalartype(x),cod,dom)

    return _is_prop_util(x, id; kwargs...)
end

isid(x::Number; kwargs...) = true, x

function isid(x::MPSBondTensor; kwargs...)
	(domain(x) == codomain(x)) || return false, 0.0
	id = isomorphism(scalartype(x), codomain(x), domain(x))
	return _is_prop_util(x, id; kwargs...)
end

function _is_prop_util(x, a; atol::Real=1.0e-14) 
	scal = dot(a,x)/dot(a,a)
	diff = x-scal*a
	scal = (scal ≈ 0.0) ? 0.0 : scal #shouldn't be necessary (and I don't think it is)
	return norm(diff)<atol,scal
end
