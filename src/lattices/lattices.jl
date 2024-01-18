include("grassmannordering.jl")


abstract type AbstractGrassmannLattice{O <: GrassmannOrdering} end

index(x::AbstractGrassmannLattice, args...; kwargs...) = error("index not implemented for grassmann lattice type $(typeof(x))")
OrderingStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = O
ConjugationStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = ConjugationStyle(O)
LayoutStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = LayoutStyle(O)
ConjugationStyle(x::AbstractGrassmannLattice) = ConjugationStyle(typeof(x))
LayoutStyle(x::AbstractGrassmannLattice) = LayoutStyle(typeof(x))

const DefaultIntegrationTruncation = truncdimcutoff(D=10000, ϵ=1.0e-12, add_back=0)

include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")
# more complicatd integration
include("integrate/integrate.jl")
include("parallelrun.jl")


function GrassmannLattice(; contour::Symbol, kwargs...)
	(contour in (:real, :imag, :mixed, :Keldysh)) || throw(ArgumentError("contour must be :real, :imag or :mixed (equivalentlt :Keldysh)"))
	if contour == :real
		return RealGrassmannLattice(; kwargs...)
	elseif contour == :imag
		return ImagGrassmannLattice(; kwargs...)
	else
		return MixedGrassmannLattice(; kwargs...)
	end
end

vacuumstate(x::AbstractGrassmannLattice) = GrassmannMPS(scalartype(x), length(x))