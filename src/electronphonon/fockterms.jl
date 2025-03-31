abstract type AbstractNTerm end

struct NTerm{N, T <: Number} <: AbstractNTerm
	positions::NTuple{N, Int}
	coeff::T

function NTerm(positions::NTuple{N, Int}; coeff::Number) where {N}
	(length(Set(positions)) == N) || throw(ArgumentError("multiple n̂ on the same position not allowed"))
	p = sortperm([positions...])
	positions = positions[p]
	new{N, typeof(coeff)}(positions, coeff)
end
end

TK.scalartype(::Type{NTerm{N, T}}) where {N, T} = T
DMRG.positions(x::NTerm) = x.positions
NTerm(a::Vararg{Int}; kwargs...) = NTerm(a; kwargs...)


function Base.convert(::Type{<:PartialDenseMPO}, x::NTerm)
	nop2 = Matrix{scalartype(x)}([0 0; 1 1])
	nop = reshape(nop2, 1,2,1,2)

	return PartialDenseMPO([nop for i in 1:length(positions(x))], [positions(x)...]) * x.coeff
end


struct ExpNTerm{N, T <:Number} <: AbstractNTerm
	data::NTerm{N, T}
end
TK.scalartype(::Type{ExpNTerm{N, T}}) where {N, T} = T
DMRG.positions(x::ExpNTerm) = x.data.positions
Base.exp(x::NTerm{N}) where N = ExpNTerm(x)
function Base.convert(::Type{<:PartialDenseMPO}, x::ExpNTerm)
	t1 = convert(PartialDenseMPO, x.data)
	return t1 + id(t1)
end


function DMRG.apply!(x::Union{NTerm, ExpNTerm}, mps::FockMPS)
    all(s -> 1 <= s <= length(mps), positions(x)) || throw(BoundsError())
    t = convert(PartialDenseMPO, x)
    apply!(t, mps)
    return mps
end
