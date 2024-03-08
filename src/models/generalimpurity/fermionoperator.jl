abstract type AbstractFTerm end

DMRG.positions(x::AbstractFTerm) = x.positions

"""
	struct TunnelingTerm <: AbstractFTerm

Fermionic twobody term
"""
struct TunnelingTerm{T <: Number} <: AbstractFTerm
	positions::Tuple{Int, Int}
	coeff::T
end

TunnelingTerm(pos::Tuple{Int, Int}; coeff::Number=1) = TunnelingTerm(pos, float(coeff))
TunnelingTerm(i::Int, j::Int; kwargs...) = TunnelingTerm((i, j); kwargs...)

tunneling(i::Int, j::Int; kwargs...) = TunnelingTerm(i, j; kwargs...)

TK.scalartype(::Type{TunnelingTerm{T}}) where {T} = T

function Base.adjoint(x::TunnelingTerm)
	i, j = positions(x)
	return TunnelingTerm((j, i), coeff=conj(x.coeff))
end

Base.copy(x::TunnelingTerm) = TunnelingTerm(positions(x), copy(x.coeff))

Base.:*(s::TunnelingTerm, m::Number) = TunnelingTerm(positions(s), coeff=s.coeff * m)
Base.:*(m::Number, s::AbstractFTerm) = s * m
Base.:/(s::AbstractFTerm, m::Number) = s * (1 / m)
Base.:+(s::AbstractFTerm) = s
Base.:-(s::AbstractFTerm) = (-1) * s


"""
	struct FourBodyTerm <: AbstractFTerm

Fermionic fourbody term
"""
struct InteractionTerm{T<:Number} <: AbstractFTerm
	positions::NTuple{4, Int}
	coeff::T
end

InteractionTerm(pos::NTuple{4, Int}; coeff::Real=1) = InteractionTerm(pos, float(coeff))
InteractionTerm(i::Int, j::Int, k::Int, l::Int; kwargs...) = InteractionTerm((i, j, k, l); kwargs...)

interaction(i::Int, j::Int, k::Int, l::Int; kwargs...) = InteractionTerm(i, j, k, l; kwargs...)

TK.scalartype(::Type{InteractionTerm{T}}) where {T} = T

function Base.adjoint(x::InteractionTerm)
	i, j, k, l = positions(x)
	return interaction((l,k,j,i), coeff=conj(x.coeff))
end

Base.:*(s::InteractionTerm, m::Number) = InteractionTerm(positions(s), coeff=s.coeff * m)


struct ImpurityHamiltonian <: AbstractImpurityHamiltonian
	data::Vector{AbstractFTerm}
	bands::Int

function ImpurityHamiltonian(data::Vector{<:AbstractFTerm}, bands::Int)
	@boundscheck begin
		for f in data
			for j in positions(f)
				(1 <= j <= bands) || throw(BoundsError(1:bands, j))
			end
		end
	end
	new(convert(Vector{AbstractFTerm}, data), bands)
end

end
ImpurityHamiltonian(data::Vector{<:AbstractFTerm}; bands::Int=1) = ImpurityHamiltonian(data, bands)
ImpurityHamiltonian(; bands::Int = 1) = ImpurityHamiltonian(Vector{AbstractFTerm}(), bands)
function Base.push!(x::ImpurityHamiltonian, f::AbstractFTerm)
	@boundscheck begin
		for j in positions(f)
			(1 <= j <= x.bands) || throw(BoundsError(1:x.bands, j))
		end
	end
	for t in x.data
		if positions(t) == positions(f)
			throw(ArgumentError("duplicate positions not allowed"))
		end
	end
	push!(x.data, f)
end 

function TK.scalartype(h::ImpurityHamiltonian)
	T = Float64
	for f in h.data
		T = promote_type(T, scalartype(f))
	end
	return T
end

