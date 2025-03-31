abstract type RealFockLattice{O<:RealFockOrdering} <: AbstractFockLattice{O} end
TK.scalartype(::Type{<:RealFockLattice}) = ComplexF64
branches(::Type{<:RealFockLattice}) = (:+, :-)
TimeOrderingStyle(x::RealFockLattice) = RealTimeOrderingStyle(x)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band
# k+1 due to the grassmann number on the boundary for the final trace

"""
	struct RealFockLattice1Order <: RealFockLattice

First order splitting of the real-time contour
"""
struct RealFockLattice1Order{O<:RealFockOrdering} <: RealFockLattice{O}
	δt::Float64
	bands::Int
	N::Int
	ordering::O

	RealFockLattice1Order(δt::Real, bands::Int, N::Int, ordering::RealFockOrdering) = new{typeof(ordering)}(convert(Float64, δt), bands, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
RealFockLattice1Order(; δt::Real, N::Int, bands::Int=1, ordering::RealFockOrdering=M1m1N1n1()) = RealFockLattice1Order(δt, bands, N, ordering)
Base.similar(x::RealFockLattice1Order; δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N, ordering::RealFockOrdering=x.ordering) = RealFockLattice1Order(δt, bands, N, ordering)
similargrassmannlattice(x::RealFockLattice1Order, δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N, 
						ordering::RealGrassmannOrdering=similargrassmannordering(x.ordering)) = GrassmannLattice(contour=:real, δt=δt, N=N, bands=bands, ordering=ordering)


function Base.getproperty(x::RealFockLattice, s::Symbol)
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
Base.length(x::RealFockLattice) = 2*x.bands * x.N

function RealFockLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return RealFockLattice1Order(; kwargs...)
	else
		error("Second orderr RealFockLattice not implemented")
	end
end


function index(x::RealFockLattice{<:M1m1N1n1}, i::Int; branch::Symbol=:+, band::Int=1)
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
function indexmappings(lattice::RealFockLattice)
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

