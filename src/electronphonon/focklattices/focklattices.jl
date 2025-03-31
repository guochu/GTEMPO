include("fockordering.jl")

abstract type AbstractFockLattice{O <: FockOrdering} end

index(x::AbstractFockLattice, args...; kwargs...) = error("index not implemented for grassmann lattice type $(typeof(x))")
OrderingStyle(::Type{<:AbstractFockLattice{O}}) where O = O()
LayoutStyle(::Type{<:AbstractFockLattice{O}}) where O = LayoutStyle(O)
ImaginaryTimeOrderingStyle(::Type{<:AbstractFockLattice{O}}) where O = ImaginaryTimeOrderingStyle(O)
RealTimeOrderingStyle(::Type{<:AbstractFockLattice{O}}) where O = RealTimeOrderingStyle(O)
OrderingStyle(x::AbstractFockLattice) = OrderingStyle(typeof(x))
ConjugationStyle(x::AbstractFockLattice) = ConjugationStyle(typeof(x))
LayoutStyle(x::AbstractFockLattice) = LayoutStyle(typeof(x))
ImaginaryTimeOrderingStyle(x::AbstractFockLattice) = ImaginaryTimeOrderingStyle(typeof(x))
RealTimeOrderingStyle(x::AbstractFockLattice) = RealTimeOrderingStyle(typeof(x))
branches(lattice::AbstractFockLattice) = branches(typeof(lattice))

include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


function FockLattice(; contour::Symbol, kwargs...)
	(contour in (:real, :imag, :Keldysh, :mixed, :KadanoffBaym)) || throw(ArgumentError("contour must be :real (equivalentlt :Keldysh), :imag or :mixed (equivalentlt :KadanoffBaym)"))
	if (contour == :real) || (contour == :Keldysh)
		return RealFockLattice(; kwargs...)
	elseif contour == :imag
		return ImagFockLattice(; kwargs...)
	else
		return MixedFockLattice(; kwargs...)
	end
end


vacuumstate(x::AbstractFockLattice) = FockMPS(scalartype(x), length(x))