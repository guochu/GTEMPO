abstract type MixedGrassmannLattice{O<:MixedGrassmannOrdering} <: AbstractGrassmannLattice{O} end
TK.scalartype(::Type{<:MixedGrassmannLattice}) = ComplexF64

struct MixedGrassmannLattice1Order{O<:MixedGrassmannOrdering} <: MixedGrassmannLattice{O}
	β::Float64
	δτ::Float64
	t::Float64
	δt::Float64
	bands::Int
	N::Int
	ordering::O
end
