abstract type AbstractRealFockLattice{O<:RealFockOrdering} <: AbstractFockLattice{O} end
TK.scalartype(::Type{<:AbstractRealFockLattice}) = ComplexF64
branches(::Type{<:AbstractRealFockLattice}) = (:+, :-)
TimeOrderingStyle(x::AbstractRealFockLattice) = RealTimeOrderingStyle(x)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band
# k+1 due to the grassmann number on the boundary for the final trace

"""
	struct RealFockLattice <: AbstractRealFockLattice

First order splitting of the real-time contour
"""
struct RealFockLattice{O<:RealFockOrdering} <: AbstractRealFockLattice{O}
	δt::Float64
	bands::Int
	N::Int
	ordering::O

	RealFockLattice(δt::Real, bands::Int, N::Int, ordering::RealFockOrdering) = new{typeof(ordering)}(convert(Float64, δt), bands, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
RealFockLattice(; δt::Real, N::Int, bands::Int=1, ordering::RealFockOrdering=M1m1N1n1()) = RealFockLattice(δt, bands, N, ordering)
Base.similar(x::RealFockLattice; δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N, ordering::RealFockOrdering=x.ordering) = RealFockLattice(δt, bands, N, ordering)
similargrassmannlattice(x::RealFockLattice, δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N,
						ordering::RealGrassmannOrdering=similargrassmannordering(x.ordering)) = GrassmannLattice(contour=:real, δt=δt, N=N, bands=bands, ordering=ordering)


function Base.getproperty(x::AbstractRealFockLattice, s::Symbol)
	if s == :t
		return x.N * x.δt
	elseif s == :ts
		return 0:x.δt:x.N*x.δt
	elseif s == :Nt
		return x.N
	else
		getfield(x, s)
	end
end
Base.length(x::AbstractRealFockLattice) = 2*x.bands * x.N


function index(x::AbstractRealFockLattice{<:M1m1N1n1}, i::Int; branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(1 <= i <= x.N) || throw(BoundsError(1:x.N, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	ifelse(branch == :+, TL-2i*bands+1+2(band-1), TL-2i*bands+2band)
end


# key is timestep, conj, branch, band
function indexmappings(lattice::AbstractRealFockLattice)
	r = Dict{Tuple{Int, Symbol, Int}, Int}()
	for i in 1:lattice.N
		for band in 1:lattice.bands
			for f in (:+, :-)
				r[(i, f, band)] = index(lattice, i, branch=f, band=band)
			end
		end
	end
	return r
end
