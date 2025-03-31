abstract type MixedFockLattice{O<:MixedFockOrdering} <: AbstractFockLattice{O} end
TK.scalartype(::Type{<:MixedFockLattice}) = ComplexF64
branches(::Type{<:MixedFockLattice}) = (:+, :-, :τ)

"""
	struct MixedGrassmannLattice1Order <: MixedGrassmannLattice

First order splitting of the real-time contour
"""
struct MixedFockLattice1Order{O<:MixedFockOrdering} <: MixedFockLattice{O}
	δt::Float64
	Nt::Int
	δτ::Float64
	Nτ::Int
	bands::Int
	ordering::O

	MixedFockLattice1Order(δt::Real, N::Int, δτ::Real, Ni::Int, bands::Int, ordering::MixedFockOrdering) = new{typeof(ordering)}(
								convert(Float64, δt), N, convert(Float64, δτ), Ni, bands, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
MixedFockLattice1Order(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, bands::Int=1, ordering::MixedFockOrdering=M1N1_m1M1n1N1m2M2n2N2()) = MixedFockLattice1Order(
							δt, Nt, δτ, Nτ, bands, ordering)
Base.similar(x::MixedFockLattice1Order; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, bands::Int=x.bands, ordering::MixedFockOrdering=x.ordering) = MixedFockLattice1Order(
			δt, Nt, δτ, Nτ, bands, ordering)
similargrassmannlattice(x::MixedFockLattice1Order, δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, bands::Int=x.bands, 
						ordering::MixedGrassmannOrdering=similargrassmannordering(x.ordering)) = GrassmannLattice(contour=:mixed, δt=δt, Nt=Nt, δτ=δτ, Nτ=Nτ, bands=bands, ordering=ordering)


function MixedFockLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return MixedFockLattice1Order(; kwargs...)
	else
		error("Second orderr MixedGrassmannLattice not implemented")
	end
end

function Base.getproperty(x::MixedFockLattice1Order, s::Symbol)
	if s == :t
		return x.Nt * x.δt
	elseif s == :β
		return x.Nτ * x.δτ
	elseif s == :T
		return 1 / x.β
	elseif s == :ts
		return 0:x.δt:x.t
	elseif s == :τs
		return 0:x.δτ:x.β
	else
		getfield(x, s)
	end
end

Base.length(x::MixedFockLattice1Order) = 2*x.bands * x.Nt + x.bands * x.Nτ

# acending order for real branch, descending order for imag time
function index(x::MixedFockLattice1Order{<:M1N1_m1M1n1N1m2M2n2N2}, i::Int; branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if branch == :τ
			(1 <= i <= x.Nτ) || throw(BoundsError(1:x.Nτ, i))
		else
			(1 <= i <= x.Nt) || throw(BoundsError(1:x.Nt, i))
		end
	end

	bands = x.bands

	if branch == :+
		2*(i-1)*bands+2+2*(band-1) + bands*x.Nτ
	elseif branch == :-
		2*(i-1)*bands+1+2*(band-1) + bands*x.Nτ
	else
		(x.Nτ-i)*bands + band
	end
end



# key is timestep, conj, branch, band
function indexmappings(lattice::MixedFockLattice1Order)
	r = Dict{Tuple{Int, Symbol, Int}, Int}()
	for i in 1:lattice.Nτ
		for band in 1:lattice.bands
			f = :τ
			r[(i, f, band)] = index(lattice, i, branch=f, band=band)
		end
	end
	for i in 1:lattice.Nt
		for band in 1:lattice.bands
			for f in (:+, :-)
				r[(i, f, band)] = index(lattice, i, branch=f, band=band)
			end
		end
	end
	return r
end

