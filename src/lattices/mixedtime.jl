abstract type MixedGrassmannLattice{O<:MixedGrassmannOrdering} <: AbstractGrassmannLattice{O} end
TK.scalartype(::Type{<:MixedGrassmannLattice}) = ComplexF64


"""
	struct MixedGrassmannLattice1Order <: MixedGrassmannLattice

First order splitting of the real-time contour
"""
struct MixedGrassmannLattice1Order{O<:MixedGrassmannOrdering} <: MixedGrassmannLattice{O}
	δt::Float64
	Nt::Int
	δτ::Float64
	Nτ::Int
	bands::Int
	ordering::O

	MixedGrassmannLattice1Order(δt::Real, N::Int, δτ::Real, Ni::Int, bands::Int, ordering::MixedGrassmannOrdering) = new{typeof(ordering)}(
								convert(Float64, δt), N, convert(Float64, δτ), Ni, bands, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
MixedGrassmannLattice1Order(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, bands::Int=1, ordering::MixedGrassmannOrdering=AABB_aaAAbbBB()) = MixedGrassmannLattice1Order(
							δt, Nt, δτ, Nτ, bands, ordering)
Base.similar(x::MixedGrassmannLattice1Order; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, bands::Int=x.bands, ordering::MixedGrassmannOrdering=x.ordering) = MixedGrassmannLattice1Order(
			δt, Nt, δτ, Nτ, bands, ordering)

function MixedGrassmannLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return MixedGrassmannLattice1Order(; kwargs...)
	else
		error("Second orderr MixedGrassmannLattice not implemented")
	end
end

function Base.getproperty(x::MixedGrassmannLattice1Order, s::Symbol)
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

Base.length(x::MixedGrassmannLattice1Order) = 4*x.bands * (x.Nt+1) + 2 * x.bands + 2*x.bands * x.Nτ

# acending order for real branch, descending order for imag time
function index(x::MixedGrassmannLattice1Order{<:A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(ArgumentError("band $band out of range"))
		(branch in (:+, :-, :τ)) || throw(BoundsError("branch must be one of :+, :- or :τ"))
		if i != 0
			if branch == :τ
				(1 <= i <= x.Nτ) || throw(BoundsError("imag time step $i out of range"))
			else
				(1 <= i <= x.Nt + 1) || throw(ArgumentError("real time step $i out of range"))
			end
		end
	end

	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1) 
	else
		k = x.Nt + 1
		if branch == :+
			ifelse(conj, 4*(i-1)*bands+2+4*(band-1), 4*(i-1)*bands+1+4*(band-1)) + 2*x.bands*(x.Nτ+1)
		elseif branch == :-
			ifelse(conj, 4*(i-1)*bands+4+4*(band-1), 4*(i-1)*bands+3+4*(band-1)) + 2*x.bands*(x.Nτ+1)
		else
			ifelse(conj, (x.Nτ-i)*2*bands + 2*band, (x.Nτ-i)*2*bands + 2*band-1) + 2*x.bands
		end
	end
end

function index(x::MixedGrassmannLattice1Order{<:A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(ArgumentError("band $band out of range"))
		(branch in (:+, :-, :τ)) || throw(BoundsError("branch must be one of :+, :- or :τ"))
		if i != 0
			if branch == :τ
				(1 <= i <= x.Nτ) || throw(BoundsError("imag time step $i out of range"))
			else
				(1 <= i <= x.Nt + 1) || throw(ArgumentError("real time step $i out of range"))
			end
		end
	end

	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1) 
	else
		k = x.Nt + 1
		if branch == :-
			ifelse(conj, 4*(i-1)*bands+2+4*(band-1), 4*(i-1)*bands+1+4*(band-1)) + 2*x.bands*(x.Nτ+1)
		elseif branch == :+
			ifelse(conj, 4*(i-1)*bands+4+4*(band-1), 4*(i-1)*bands+3+4*(band-1)) + 2*x.bands*(x.Nτ+1)
		else
			ifelse(conj, (x.Nτ-i)*2*bands + 2*band, (x.Nτ-i)*2*bands + 2*band-1) + 2*x.bands
		end
	end
end

# acending order for real branch, descending order for imag time
function index(x::MixedGrassmannLattice1Order{<:A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(ArgumentError("band $band out of range"))
		(branch in (:+, :-, :τ)) || throw(BoundsError("branch must be one of :+, :- or :τ"))
		if i != 0
			if branch == :τ
				(1 <= i <= x.Nτ) || throw(BoundsError("imag time step $i out of range"))
			else
				(1 <= i <= x.Nt + 1) || throw(ArgumentError("real time step $i out of range"))
			end
		end
	end
	
	bands = x.bands
	if i == 0
		ifelse(conj, 2*bands+1-band, band) 
	else
		k = x.Nt + 1
		if branch == :+
			ifelse(conj, 2*bands*(k-i) + 2*bands-band+1, 2*bands*(k-i) + band ) + 2*x.bands*(x.Nτ+1)
		elseif branch == :-
			ifelse(conj, 2*bands*(i-1) + 2*bands-band+1, 2*bands*(i-1) + band ) + 2 * bands * k + 2*x.bands*(x.Nτ+1)
		else
			ifelse(conj, (x.Nτ-i)*2*bands + 2bands+1-band, (x.Nτ-i)*2*bands + band) + 2*x.bands
		end
	end
end


# key is timestep, conj, branch, band
function indexmappings(lattice::MixedGrassmannLattice1Order)
	r = Dict{Tuple{Int, Bool, Symbol, Int}, Int}()
	for i in 1:lattice.Nτ
		for c in (true, false)
			for band in 1:lattice.bands
				f = :τ
				r[(i, c, f, band)] = index(lattice, i, conj=c, branch=f, band=band)
			end
		end
	end
	for i in 0:lattice.Nt+1
		for c in (true, false)
			for band in 1:lattice.bands
				if i == 0
					r[(i, c, :+, band)] = index(lattice, i, conj=c, branch=:+, band=band)
				else
					for f in (:+, :-)
						r[(i, c, f, band)] = index(lattice, i, conj=c, branch=f, band=band)
					end
				end
			end
		end
	end
	return r
end

