abstract type ImagGrassmannLattice{O<:ImagGrassmannOrdering} <: AbstractGrassmannLattice{O} end
TK.scalartype(::Type{<:ImagGrassmannLattice}) = Float64
branches(::Type{<:ImagGrassmannLattice}) = (:τ,)
TimeOrderingStyle(x::ImagGrassmannLattice) = ImaginaryTimeOrderingStyle(x)

struct ImagGrassmannLattice1Order{O<:ImagGrassmannOrdering} <: ImagGrassmannLattice{O}
	δτ::Float64
	bands::Int
	N::Int
	ordering::O

	ImagGrassmannLattice1Order(δτ::Real, bands::Int, N::Int, ordering::ImagGrassmannOrdering) = new{typeof(ordering)}(convert(Float64, δτ), bands, N, ordering)
end

ImagGrassmannLattice1Order(; δτ::Real, N::Int, bands::Int=1, ordering::ImagGrassmannOrdering=A1Ā1B1B̄1()) = ImagGrassmannLattice1Order(δτ, bands, N, ordering)
Base.similar(x::ImagGrassmannLattice1Order; δτ::Real=x.δτ, bands::Int=x.bands, N::Int=x.N, ordering::ImagGrassmannOrdering=x.ordering) = ImagGrassmannLattice1Order(δτ, bands, N, ordering)

function ImagGrassmannLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return ImagGrassmannLattice1Order(; kwargs...)
	else
		error("Second orderr ImagGrassmannLattice not implemented")
	end
end

Base.length(x::ImagGrassmannLattice) = 2*x.bands * (x.k+1)
function Base.getproperty(x::ImagGrassmannLattice1Order, s::Symbol)
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
	else
		getfield(x, s)
	end
end
# function Base.:(==)(a::ImagGrassmannLattice1Order, b::ImagGrassmannLattice1Order)
# 	return (a.δτ == b.δτ) && (a.bands == b.bands) && (a.N == b.N) && (OrderingStype(a) == OrderingStyle(b))
# end


# ab\bar{b}\bar{a} a_2b_2\bar{b}_2\bar{a}_2 a_1b_1\bar{b}_1\bar{a}_1
function index(x::ImagGrassmannLattice{<:A1B1B̄1Ā1}, i::Int; conj::Bool, band::Int=1, branch::Symbol=:τ)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	# TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2bands+1-band, band)
	else
		ifelse(conj, (x.k-i)*2*bands + 2bands+1-band, (x.k-i)*2*bands + band) + 2*bands
	end	
end
# a\bar{a}b\bar{b} a_2\bar{a}_2b_2\bar{b}_2 a_1\bar{a}_1b_1\bar{b}_1
function index(x::ImagGrassmannLattice{<:A1Ā1B1B̄1}, i::Int; conj::Bool, band::Int=1, branch::Symbol=:τ)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	# TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1)
	else
		ifelse(conj, (x.k-i)*2*bands + 2*band, (x.k-i)*2*bands + 2*band-1) + 2*bands
	end	
end
# a\bar{a}b\bar{b} a_2\bar{a}_2a_1\bar{a}_1 b_2\bar{b}_2b_1\bar{b}_1
function index(x::ImagGrassmannLattice{<:A2Ā2A1Ā1B2B̄2B1B̄1}, i::Int; conj::Bool, band::Int=1, branch::Symbol=:τ)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	n = 2 * x.k 
	if i == 0
		ifelse(conj, 2*band, 2*band-1)
	else
		ifelse(conj, (band-1) * n + (x.k-i)*2 + 2, (band-1) * n + (x.k-i)*2 + 1) + 2 * x.bands
	end	
end

function index(x::ImagGrassmannLattice{<:Ā2A1B̄2B1}, i::Int; conj::Bool, band::Int=1, branch::Symbol=:τ)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	n = 2 * x.k 
	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1)
	else
		if (i == 1) && conj
			length(x) - bands + band
		elseif (i == x.k) && (!conj)
			2*bands + band
		else
			ifelse(conj, (x.k-i)*2*bands + 2*(band-1)+1, (x.k-i-1)*2*bands + 2*band) + 3*bands
		end
	end	
end


function indexmappings(lattice::ImagGrassmannLattice)
	r = Dict{Tuple{Int, Bool, Symbol, Int}, Int}()
	for i in 0:lattice.k
		for c in (true, false)
			for band in 1:lattice.bands
				r[(i, c, :τ, band)] = index(lattice, i, conj=c, band=band)
			end
		end
	end
	return r
end

function band_boundary(lattice::ImagGrassmannLattice{<:A1B1B̄1Ā1}, j::Int)
    posa = index(lattice, j, conj=false, band=1)
    posb = index(lattice, j, conj=true, band=1)     
    return posa, posb
end
function band_boundary(lattice::ImagGrassmannLattice{<:A1Ā1B1B̄1}, j::Int)
    posa = index(lattice, j, conj=false, band=1)
    posb = index(lattice, j, conj=true, band=lattice.bands)     
    return posa, posb
end
