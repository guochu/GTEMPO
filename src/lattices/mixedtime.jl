abstract type AbstractMixedGrassmannLattice{O<:MixedGrassmannOrdering} <: AbstractGrassmannLattice{O} end
TK.scalartype(::Type{<:AbstractMixedGrassmannLattice}) = ComplexF64
branches(::Type{<:AbstractMixedGrassmannLattice}) = (:+, :-, :τ)

"""
	struct MixedGrassmannLattice <: AbstractMixedGrassmannLattice

Mixed (Kadanoff–Baym) Grassmann lattice.
"""
struct MixedGrassmannLattice{O<:MixedGrassmannOrdering} <: AbstractMixedGrassmannLattice{O}
	δt::Float64
	Nt::Int
	δτ::Float64
	Nτ::Int
	bands::Int
	ordering::O

	MixedGrassmannLattice(δt::Real, N::Int, δτ::Real, Ni::Int, bands::Int, ordering::MixedGrassmannOrdering) = new{typeof(ordering)}(
								convert(Float64, δt), N, convert(Float64, δτ), Ni, bands, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
MixedGrassmannLattice(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, bands::Int=1, ordering::MixedGrassmannOrdering=A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2()) = MixedGrassmannLattice(
							δt, Nt, δτ, Nτ, bands, ordering)
Base.similar(x::MixedGrassmannLattice; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ, bands::Int=x.bands, ordering::MixedGrassmannOrdering=x.ordering) = MixedGrassmannLattice(
			δt, Nt, δτ, Nτ, bands, ordering)

function Base.getproperty(x::MixedGrassmannLattice, s::Symbol)
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

Base.length(x::MixedGrassmannLattice) = 4*x.bands * x.kt + 2*x.bands * x.kτ

# acending order for real branch, descending order for imag time
function _index(x::MixedGrassmannLattice{<:A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
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
	if branch == :+
		ifelse(conj, 4*(i-1)*bands+2+4*(band-1), 4*(i-1)*bands+1+4*(band-1)) + 2*x.bands*x.kτ
	elseif branch == :-
		ifelse(conj, 4*(i-1)*bands+4+4*(band-1), 4*(i-1)*bands+3+4*(band-1)) + 2*x.bands*x.kτ
	else
		ifelse(conj, (x.Nτ+1-i)*2*bands + 2*band, (x.Nτ+1-i)*2*bands + 2*band-1)
	end
end

function _index(x::MixedGrassmannLattice{<:A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
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
	if branch == :-
		ifelse(conj, 4*(i-1)*bands+2+4*(band-1), 4*(i-1)*bands+1+4*(band-1)) + 2*x.bands*x.kτ
	elseif branch == :+
		ifelse(conj, 4*(i-1)*bands+4+4*(band-1), 4*(i-1)*bands+3+4*(band-1)) + 2*x.bands*x.kτ
	else
		ifelse(conj, (x.Nτ+1-i)*2*bands + 2*band, (x.Nτ+1-i)*2*bands + 2*band-1)
	end
end


# internal site-position map; the public accessor is
# `lattice[GrassmannContourIndex(i, branch, band, conj)]`
function _index(x::MixedGrassmannLattice{<:Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
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
	if branch == :τ
		if (i == 1) && conj
			bands + band
		elseif (i == x.kτ) && (!conj)
			band
		else
			ifelse(conj, (x.kτ-i)*2*bands + 2*(band-1)+1, (x.kτ-i-1)*2*bands + 2*band) + 4bands
		end
	else
		n_imag = 2bands*x.kτ
		kt = x.kt
		if (i == 1) && conj && (branch==:+)
			2*band + 2bands
		elseif (i == 1) && (!conj) && (branch==:-)
			2*(band-1)+1 + 2bands
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

# A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1: imaginary part A1B1B̄1Ā1, real part
# interleaves the two branches (canonical ordering for the fast propagator)
function _index(x::MixedGrassmannLattice{<:A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
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
	TL = length(x)
	bands = x.bands
	if branch == :-
		ifelse(conj, 4i*bands+2*bands+band, 4i*bands+band) + (x.Nτ-1)*2*bands
	elseif branch == :+
		ifelse(conj, 4i*bands+bands+band, 4i*bands+3*bands+band) + (x.Nτ-1)*2*bands
	else
		ifelse(conj, (x.Nτ+1-i)*2*bands + 2bands+1-band, (x.Nτ+1-i)*2*bands + band)
	end
end

# key is timestep, conj, branch, band
function indexmappings(lattice::MixedGrassmannLattice)
	r = Dict{Tuple{Int, Bool, Symbol, Int}, Int}()
	for i in 1:lattice.kτ
		for c in (true, false)
			for band in 1:lattice.bands
				f = :τ
				r[(i, c, f, band)] = lattice[GrassmannContourIndex(i, f, band, c)]
			end
		end
	end
	for i in 1:lattice.kt
		for c in (true, false)
			for band in 1:lattice.bands
				for f in (:+, :-)
					r[(i, c, f, band)] = lattice[GrassmannContourIndex(i, f, band, c)]
				end
			end
		end
	end
	return r
end

