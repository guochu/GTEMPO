abstract type AbstractImagGrassmannLattice{O<:ImagGrassmannOrdering} <: AbstractGrassmannLattice{O} end
TK.scalartype(::Type{<:AbstractImagGrassmannLattice}) = Float64
branches(::Type{<:AbstractImagGrassmannLattice}) = (:τ,)
TimeOrderingStyle(x::AbstractImagGrassmannLattice) = ImaginaryTimeOrderingStyle(x)

struct ImagGrassmannLattice{O<:ImagGrassmannOrdering} <: AbstractImagGrassmannLattice{O}
	δτ::Float64
	bands::Int
	N::Int
	ordering::O

	ImagGrassmannLattice(δτ::Real, bands::Int, N::Int, ordering::ImagGrassmannOrdering) = new{typeof(ordering)}(convert(Float64, δτ), bands, N, ordering)
end

ImagGrassmannLattice(; δτ::Real, N::Int, bands::Int=1, ordering::ImagGrassmannOrdering=A1Ā1B1B̄1()) = ImagGrassmannLattice(δτ, bands, N, ordering)
Base.similar(x::ImagGrassmannLattice; δτ::Real=x.δτ, bands::Int=x.bands, N::Int=x.N, ordering::ImagGrassmannOrdering=x.ordering) = ImagGrassmannLattice(δτ, bands, N, ordering)

Base.length(x::AbstractImagGrassmannLattice) = 2*x.bands * x.k
function Base.getproperty(x::ImagGrassmannLattice, s::Symbol)
	if (s == :k) || ( s == :kτ)
		return x.N+1
	elseif s == :τs
		return 0:x.δτ:x.N*x.δτ
	elseif s == :β
		return x.N * x.δτ
	elseif s == :T
		return 1 / x.β
	elseif s == :Nτ
		return x.N
	elseif s == :orbitals
		bands = x.bands
		iseven(bands) || throw(ArgumentError("even number of bands expected"))
		return div(bands, 2)
	else
		getfield(x, s)
	end
end
# function Base.:(==)(a::ImagGrassmannLattice, b::ImagGrassmannLattice)
# 	return (a.δτ == b.δτ) && (a.bands == b.bands) && (a.N == b.N) && (OrderingStype(a) == OrderingStyle(b))
# end


# ab\bar{b}\bar{a} a_2b_2\bar{b}_2\bar{a}_2 a_1b_1\bar{b}_1\bar{a}_1
function _index(x::AbstractImagGrassmannLattice{<:A1B1B̄1Ā1}, i::Int; conj::Bool, band::Int=1, branch::Symbol=:τ)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	# TL = length(x)
	bands = x.bands
	ifelse(conj, (x.k-i)*2*bands + 2bands+1-band, (x.k-i)*2*bands + band)
end
# a\bar{a}b\bar{b} a_2\bar{a}_2b_2\bar{b}_2 a_1\bar{a}_1b_1\bar{b}_1
function _index(x::AbstractImagGrassmannLattice{<:A1Ā1B1B̄1}, i::Int; conj::Bool, band::Int=1, branch::Symbol=:τ)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	# TL = length(x)
	bands = x.bands
	ifelse(conj, (x.k-i)*2*bands + 2*band, (x.k-i)*2*bands + 2*band-1)
end

function _index(x::AbstractImagGrassmannLattice{<:Ā2A1B̄2B1}, i::Int; conj::Bool, band::Int=1, branch::Symbol=:τ)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	bands = x.bands
	if (i == 1) && conj
		length(x) - bands + band
	elseif (i == x.k) && (!conj)
		band
	else
		ifelse(conj, (x.k-i)*2*bands + 2*(band-1)+1, (x.k-i-1)*2*bands + 2*band) + bands
	end
end


function indexmappings(lattice::AbstractImagGrassmannLattice)
	r = Dict{Tuple{Int, Bool, Symbol, Int}, Int}()
	for i in 1:lattice.k
		for c in (true, false)
			for band in 1:lattice.bands
				r[(i, c, :τ, band)] = lattice[GrassmannContourIndex(i, :τ, band, c)]
			end
		end
	end
	return r
end

function band_boundary(lattice::AbstractImagGrassmannLattice{<:A1B1B̄1Ā1}, j::Int)
    posa = lattice[GrassmannContourIndex(j, :τ, 1, false)]
    posb = lattice[GrassmannContourIndex(j, :τ, 1, true)]
    return posa, posb
end
function band_boundary(lattice::AbstractImagGrassmannLattice{<:A1Ā1B1B̄1}, j::Int)
    posa = lattice[GrassmannContourIndex(j, :τ, 1, false)]
    posb = lattice[GrassmannContourIndex(j, :τ, lattice.bands, true)]
    return posa, posb
end
