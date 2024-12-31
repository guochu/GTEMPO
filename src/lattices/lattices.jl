include("grassmannordering.jl")
include("branch.jl")

abstract type AbstractGrassmannLattice{O <: GrassmannOrdering} end

index(x::AbstractGrassmannLattice, args...; kwargs...) = error("index not implemented for grassmann lattice type $(typeof(x))")
Base.getindex(x::AbstractGrassmannLattice, a::ContourIndex) = index(x, a.j, conj=a.conj, branch=branch(a), band=a.band)
OrderingStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = O()
ConjugationStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = ConjugationStyle(O)
LayoutStyle(::Type{<:AbstractGrassmannLattice{O}}) where O = LayoutStyle(O)
OrderingStyle(x::AbstractGrassmannLattice) = OrderingStyle(typeof(x))
ConjugationStyle(x::AbstractGrassmannLattice) = ConjugationStyle(typeof(x))
LayoutStyle(x::AbstractGrassmannLattice) = LayoutStyle(typeof(x))


include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")
# more complicatd integration
include("integrate/integrate.jl")
include("parallelrun.jl")


function GrassmannLattice(; contour::Symbol, kwargs...)
	(contour in (:real, :imag, :Keldysh, :mixed, :Kadanoff)) || throw(ArgumentError("contour must be :real (equivalentlt :Keldysh), :imag or :mixed (equivalentlt :Kadanoff)"))
	if (contour == :real) || (contour == :Keldysh)
		return RealGrassmannLattice(; kwargs...)
	elseif contour == :imag
		return ImagGrassmannLattice(; kwargs...)
	else
		return MixedGrassmannLattice(; kwargs...)
	end
end

vacuumstate(x::AbstractGrassmannLattice) = GrassmannMPS(scalartype(x), length(x))


"""
	matchindices(tsc::AbstractGrassmannLattice, src::AbstractGrassmannLattice)

Find the permutations that match src to tsc
"""
function matchindices(tsc::AbstractGrassmannLattice, src::AbstractGrassmannLattice) 
	((length(src) == length(tsc)) && (src.bands == tsc.bands)) || throw(ArgumentError("lattice size mismatch"))
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

branches(lattice::AbstractGrassmannLattice) = branches(typeof(lattice))

# This function is used in some special cases
band_boundary(lattice::AbstractGrassmannLattice, j::Int; kwargs...) = error("indexmappings not implemented for lattice type $(typeof(lattice))")

"""
	swapband(mps::GrassmannMPS, x::AbstractGrassmannLattice, b1::Int, b2::Int; kwargs...)

Swap the two bands b1 b2 of a mps
"""
swapband(mps::GrassmannMPS, x::AbstractGrassmannLattice, b1::Int, b2::Int; kwargs...) = swapband!(copy(mps), x, b1, b2; kwargs...)
swapband!(mps::GrassmannMPS, x::AbstractGrassmannLattice, b1::Int, b2::Int; kwargs...) = _swapband!(mps, x, b1, b2; kwargs...)

function _swapband!(mps, x::AbstractGrassmannLattice, b1::Int, b2::Int; kwargs...)
	(length(mps) == length(x)) || throw(DimensionMismatch("size mismatch"))
	(b1 == b2) && return mps
	perm = swapbandperm(x, b1, b2)
	return permute!(mps, [perm[i] for i in 1:length(perm)]; kwargs...)
end


function swapbandperm(x::AbstractGrassmannLattice, b1::Int, b2::Int)
	@boundscheck begin
		(1 <= b1 <= x.bands) || throw(BoundsError(1:x.bands, b1))
		(1 <= b2 <= x.bands) || throw(BoundsError(1:x.bands, b2))
	end
	change_band(b) = ifelse(b == b1, b2, ifelse(b == b2, b1, b))
	r1 = indexmappings(x)
	r2 = Dict((j, c, b, change_band(band))=>pos for ((j, c, b, band), pos) in r1)
	return Dict(r2[k1]=>v1 for (k1, v1) in r1)	
end

"""
	fillband(lattice::AbstractGrassmannLattice, gmps::GrassmannMPS; band::Int=1)
"""
function fillband(lattice::AbstractGrassmannLattice, gmps::GrassmannMPS; band::Int=1)
	data = _fillband(lattice, gmps, band=band)
	return GrassmannMPS(data, scaling=scaling(gmps)^(length(gmps)/length(lattice)))
end

function _fillband(lattice::AbstractGrassmannLattice, gmps; band::Int=1)
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	lattice2 = similar(lattice, bands=1)
	(length(lattice2) == length(gmps)) || throw(ArgumentError("GMPS size mismatch with lattice $(lattice2)"))
	r2 = indexmappings(lattice2)
	r1 = indexmappings(lattice)
	mm = Dict(r1[(j, c, b, band)]=>pos for ((j, c, b, bj), pos) in r2)

	data = similar(gmps.data, length(lattice))

	i = 1
	while isnothing(get(mm, i, nothing))
		i += 1
	end
	leftspace = space_l(gmps[mm[i]])

	for i in 1:length(lattice)
		pos2 = get(mm, i, nothing)
		if isnothing(pos2)
			data[i] = trivial_sitetenor(scalartype(gmps), leftspace)
		else
			data[i] = gmps[pos2]
			leftspace = space_r(data[i])'
		end
	end
	return data
end


function trivial_sitetenor(::Type{T}, leftspace) where {T <: Number}
	v = TensorMap(ds->zeros(T, ds), leftspace ⊗ grassmannpspace() ← leftspace )
	for s in sectors(leftspace)
		d = dim(leftspace, s)
		copy!(v[(s, Irrep[ℤ₂](0), s)], reshape(one(zeros(d, d)), d, 1, d))
	end
	return v
end


include("integrateband.jl")

