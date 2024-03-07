abstract type AbstractFTerm end

DMRG.positions(x::AbstractFTerm) = x.positions

"""
	struct TwoBodyFTerm <: AbstractFTerm

Fermionic twobody term
"""
struct TwoBodyFTerm{T <: Number} <: AbstractFTerm
	positions::Tuple{Int, Int}
	coeff::T
end

TwoBodyFTerm(pos::Tuple{Int, Int}; coeff::Number=1) = TwoBodyFTerm(pos, float(coeff))
TwoBodyFTerm(i::Int, j::Int; kwargs...) = TwoBodyFTerm((i, j); kwargs...)

twobody(i::Int, j::Int; kwargs...) = TwoBodyFTerm(i, j; kwargs...)

TK.scalartype(::Type{TwoBodyFTerm{T}}) where {T} = T

function Base.adjoint(x::TwoBodyFTerm)
	i, j = positions(x)
	return TwoBodyFTerm((j, i), coeff=conj(x.coeff))
end

Base.copy(x::TwoBodyFTerm) = TwoBodyFTerm(positions(x), copy(x.coeff))

Base.:*(s::TwoBodyFTerm, m::Number) = TwoBodyFTerm(positions(s), coeff=s.coeff * m)
Base.:*(m::Number, s::AbstractFTerm) = s * m
Base.:/(s::AbstractFTerm, m::Number) = s * (1 / m)
Base.:+(s::AbstractFTerm) = s
Base.:-(s::AbstractFTerm) = (-1) * s


"""
	struct FourBodyTerm <: AbstractFTerm

Fermionic fourbody term
"""
struct FourBodyFTerm{T<:Number} <: AbstractFTerm
	positions::NTuple{4, Int}
	coeff::T
end

FourBodyFTerm(pos::NTuple{4, Int}; coeff::Real=1) = FourBodyFTerm(pos, float(coeff))
FourBodyFTerm(i::Int, j::Int, k::Int, l::Int; kwargs...) = FourBodyFTerm((i, j, k, l); kwargs...)

fourbody(i::Int, j::Int, k::Int, l::Int; kwargs...) = FourBodyFTerm(i, j, k, l; kwargs...)

TK.scalartype(::Type{FourBodyFTerm{T}}) where {T} = T

function Base.adjoint(x::FourBodyFTerm)
	i, j, k, l = positions(x)
	return fourbody((l,k,j,i), coeff=conj(x.coeff))
end

Base.:*(s::FourBodyFTerm, m::Number) = FourBodyFTerm(positions(s), coeff=s.coeff * m)


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

