include("grassmannordering.jl")


abstract type AbstractGrassmannLattice{O <: GrassmannOrdering} end

index(x::AbstractGrassmannLattice, args...; kwargs...) = error("index not implemented for grassmann lattice type $(typeof(x))")
OrderingStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = O()
ConjugationStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = ConjugationStyle(O)
LayoutStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = LayoutStyle(O)
OrderingStyle(x::AbstractGrassmannLattice) = OrderingStyle(typeof(x))
ConjugationStyle(x::AbstractGrassmannLattice) = ConjugationStyle(typeof(x))
LayoutStyle(x::AbstractGrassmannLattice) = LayoutStyle(typeof(x))


include("imaginarytime.jl")
include("realtime.jl")
# more complicatd integration
include("integrate/integrate.jl")
include("parallelrun.jl")


function GrassmannLattice(; contour::Symbol, kwargs...)
	(contour in (:real, :imag, :Keldysh)) || throw(ArgumentError("contour must be :real (equivalentlt :Keldysh) or :imag"))
	if (contour == :real) || (contour == :Keldysh)
		return RealGrassmannLattice(; kwargs...)
	else
		return ImagGrassmannLattice(; kwargs...)
	end
end

vacuumstate(x::AbstractGrassmannLattice) = GrassmannMPS(scalartype(x), length(x))


"""
	matchindices(tsc::AbstractGrassmannLattice, src::AbstractGrassmannLattice)

Find the permutations that match src to tsc
"""
function matchindices(tsc::AbstractGrassmannLattice, src::AbstractGrassmannLattice) 
	((src.N == tsc.N) && (src.bands == tsc.bands)) || throw(ArgumentError("lattice size mismatch"))
	r1 = indexmappings(src)
	r2 = indexmappings(tsc)
	return Dict(r2[k1]=>v1 for (k1, v1) in r1)
end
function matchindices2(tsc::AbstractGrassmannLattice, src::AbstractGrassmannLattice) 
	mapping = matchindices(tsc, src)
	return [mapping[i] for i in 1:length(mapping)]
end

"""
	indexmappings(lattice::AbstractGrassmannLattice)

Return all the index mappings of a lattice
This is an internal function used for changing lattice ordering
"""
indexmappings(lattice::AbstractGrassmannLattice) = error("indexmappings not implemented for lattice type $(typeof(lattice))")


# This function is used in some special cases
band_boundary(lattice::AbstractGrassmannLattice, j::Int; kwargs...) = error("indexmappings not implemented for lattice type $(typeof(lattice))")