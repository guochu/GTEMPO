abstract type MixedGrassmannLattice{O<:MixedGrassmannOrdering} <: AbstractGrassmannLattice{O} end
TK.scalartype(::Type{<:MixedGrassmannLattice}) = ComplexF64
branches(::Type{<:MixedGrassmannLattice}) = (:+, :-, :τ)

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
MixedGrassmannLattice1Order(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, bands::Int=1, ordering::MixedGrassmannOrdering=AĀBB̄_aāAĀbb̄BB̄()) = MixedGrassmannLattice1Order(
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
	elseif s == :kt
		return x.Nt + 1
	elseif s == :kτ
		return x.Nτ + 1
	elseif s == :ts
		return 0:x.δt:x.t
	elseif s == :τs
		return 0:x.δτ:x.β
	elseif s == :orbitals
		bands = x.bands
		iseven(bands) || throw(ArgumentError("even number of bands expected"))
		return div(bands, 2)
	else
		getfield(x, s)
	end
end

Base.length(x::MixedGrassmannLattice1Order) = 4*x.bands * x.kt + 2 * x.bands + 2*x.bands * (x.Nτ+1)

# acending order for real branch, descending order for imag time
function index(x::MixedGrassmannLattice1Order{<:A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if i != 0
			if branch == :τ
				(1 <= i <= x.Nτ + 1) || throw(BoundsError(1:x.kτ, i))
			else
				(1 <= i <= x.Nt + 1) || throw(BoundsError(1:x.kt, i))
			end
		end
	end

	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1) 
	else
		k = x.Nt + 1
		if branch == :+
			ifelse(conj, 4*(i-1)*bands+2+4*(band-1), 4*(i-1)*bands+1+4*(band-1)) + 2*x.bands*(x.kτ+1)
		elseif branch == :-
			ifelse(conj, 4*(i-1)*bands+4+4*(band-1), 4*(i-1)*bands+3+4*(band-1)) + 2*x.bands*(x.kτ+1)
		else
			# if i == 1
			# 	index(x, 1, conj=conj, branch=:-, band=band)
			# else
			# 	ifelse(conj, (x.Nτ+1-i)*2*bands + 2*band, (x.Nτ+1-i)*2*bands + 2*band-1) + 2*x.bands
			# end
			ifelse(conj, (x.Nτ+1-i)*2*bands + 2*band, (x.Nτ+1-i)*2*bands + 2*band-1) + 2*x.bands
		end
	end
end

function index(x::MixedGrassmannLattice1Order{<:A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if i != 0
			if branch == :τ
				(1 <= i <= x.Nτ + 1) || throw(BoundsError(1:x.kτ, i))
			else
				(1 <= i <= x.Nt + 1) || throw(BoundsError(1:x.kt, i))
			end
		end
	end

	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1) 
	else
		k = x.Nt + 1
		if branch == :-
			ifelse(conj, 4*(i-1)*bands+2+4*(band-1), 4*(i-1)*bands+1+4*(band-1)) + 2*x.bands*(x.kτ+1)
		elseif branch == :+
			ifelse(conj, 4*(i-1)*bands+4+4*(band-1), 4*(i-1)*bands+3+4*(band-1)) + 2*x.bands*(x.kτ+1)
		else
			# if i == 1
			# 	index(x, 1, conj=conj, branch=:-, band=band)
			# else
			# 	ifelse(conj, (x.Nτ+1-i)*2*bands + 2*band, (x.Nτ+1-i)*2*bands + 2*band-1) + 2*x.bands
			# end
			ifelse(conj, (x.Nτ+1-i)*2*bands + 2*band, (x.Nτ+1-i)*2*bands + 2*band-1) + 2*x.bands
		end
	end
end

# acending order for real branch, descending order for imag time
function index(x::MixedGrassmannLattice1Order{<:A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if i != 0
			if branch == :τ
				(1 <= i <= x.Nτ + 1) || throw(BoundsError(1:x.kτ, i))
			else
				(1 <= i <= x.Nt + 1) || throw(BoundsError(1:x.kt, i))
			end
		end
	end
	
	bands = x.bands
	if i == 0
		ifelse(conj, 2*bands+1-band, band) 
	else
		k = x.Nt + 1
		if branch == :+
			ifelse(conj, 2*bands*(k-i) + 2*bands-band+1, 2*bands*(k-i) + band ) + 2*x.bands*(x.kτ+1)
		elseif branch == :-
			ifelse(conj, 2*bands*(i-1) + 2*bands-band+1, 2*bands*(i-1) + band ) + 2 * bands * k + 2*x.bands*(x.kτ+1)
		else
			# if i == 1
			# 	index(x, 1, conj=conj, branch=:-, band=band)
			# else
			# 	ifelse(conj, (x.Nτ+1-i)*2*bands + 2bands+1-band, (x.Nτ+1-i)*2*bands + band) + 2*x.bands
			# end
			ifelse(conj, (x.Nτ+1-i)*2*bands + 2bands+1-band, (x.Nτ+1-i)*2*bands + band) + 2*x.bands
		end
	end
end

"""
	index(x::MixedGrassmannLattice1Order, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
"""
function index(x::MixedGrassmannLattice1Order{<:Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if i != 0
			if branch == :τ
				(1 <= i <= x.Nτ + 1) || throw(BoundsError(1:x.kτ, i))
			else
				(1 <= i <= x.Nt + 1) || throw(BoundsError(1:x.kt, i))
			end
		end
	end
	
	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1) 
	else
		if branch == :τ
			if (i == 1) && conj
				3bands + band
			elseif (i == x.kτ) && (!conj)
				2bands + band
			else
				ifelse(conj, (x.kτ-i)*2*bands + 2*(band-1)+1, (x.kτ-i-1)*2*bands + 2*band) + 6bands
			end
		else
			n_imag = 2bands*(x.kτ+1)
			kt = x.kt			
			if (i == 1) && conj && (branch==:+)
				2*band + 4bands
			elseif (i == 1) && (!conj) && (branch==:-)
				2*(band-1)+1 + 4bands
			elseif (i == kt) && (!conj) && (branch==:+)
				length(x) - 2bands + 2*band
			elseif (i == kt) && conj && (branch==:-)
				length(x) - 2bands + 2*(band-1)+1
			else
				if branch == :-
					ifelse(conj, (i-1)*4*bands + 4*(band-1)+1, (i-2)*4*bands + 4*(band-1)+2) + n_imag + 2bands
				else
					ifelse(conj, (i-2)*4*bands + 4*(band-1)+3, (i-1)*4*bands + 4*band) + n_imag + 2bands
				end
			end
		end
	end
end

# key is timestep, conj, branch, band
function indexmappings(lattice::MixedGrassmannLattice1Order)
	r = Dict{Tuple{Int, Bool, Symbol, Int}, Int}()
	for i in 1:lattice.kτ
		for c in (true, false)
			for band in 1:lattice.bands
				f = :τ
				r[(i, c, f, band)] = index(lattice, i, conj=c, branch=f, band=band)
			end
		end
	end
	for i in 0:lattice.kt
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

