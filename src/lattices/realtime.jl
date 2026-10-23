abstract type AbstractRealGrassmannLattice{O<:RealGrassmannOrdering} <: AbstractGrassmannLattice{O} end
TK.scalartype(::Type{<:AbstractRealGrassmannLattice}) = ComplexF64
branches(::Type{<:AbstractRealGrassmannLattice}) = (:+, :-)
TimeOrderingStyle(x::AbstractRealGrassmannLattice) = RealTimeOrderingStyle(x)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band

"""
	struct RealGrassmannLattice <: AbstractRealGrassmannLattice

Real-time (Keldysh) Grassmann lattice.
"""
struct RealGrassmannLattice{O<:RealGrassmannOrdering} <: AbstractRealGrassmannLattice{O}
	δt::Float64
	bands::Int
	N::Int
	ordering::O

	RealGrassmannLattice(δt::Real, bands::Int, N::Int, ordering::RealGrassmannOrdering) = new{typeof(ordering)}(convert(Float64, δt), bands, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
RealGrassmannLattice(; δt::Real, N::Int, bands::Int=1, ordering::RealGrassmannOrdering=A1Ā1a1ā1B1B̄1b1b̄1()) = RealGrassmannLattice(δt, bands, N, ordering)
Base.similar(x::RealGrassmannLattice; δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N, ordering::RealGrassmannOrdering=x.ordering) = RealGrassmannLattice(δt, bands, N, ordering)

function Base.getproperty(x::AbstractRealGrassmannLattice, s::Symbol)
	if (s == :k) || (s == :kt)
		return x.N+1
	elseif s == :t
		return x.N * x.δt
	elseif s == :ts
		return 0:x.δt:x.N*x.δt
	elseif s == :Nt
		return x.N
	elseif s == :orbitals
		bands = x.bands
		iseven(bands) || throw(ArgumentError("even number of bands expected"))
		return div(bands, 2)
	else
		getfield(x, s)
	end
end
Base.length(x::AbstractRealGrassmannLattice) = 4*x.bands * x.k
# Base.length(x::AbstractRealGrassmannLattice) = _length(x.bands, x.k)
# _length(bands::Int, k::Int) = 4*bands*k
makestep(x::RealGrassmannLattice) = RealGrassmannLattice(x.δt, x.bands, x.N+1, x.ordering)


"""
	makestep(x::GrassmannMPS, lattice::AbstractRealGrassmannLattice)
Increase the time step of MPS-IF and the lattice by 1
"""
function makestep(lattice::AbstractRealGrassmannLattice, x::GrassmannMPS...)
	(all(v->length(v)==length(lattice), x)) || throw(DimensionMismatch())
	(LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("makestep only implemented for TimeLocalLayout"))
	# (length(x) < length(lattice)) || throw(ArgumentError("MPS-IF has reached amximum length"))
	lattice2 = makestep(lattice)
	x2 = map(y->_makestep_util(lattice2, y), x)
	return lattice2, x2...
end

function _makestep_util(lattice2::AbstractRealGrassmannLattice, x::GrassmannMPS)
	L = length(x)
	x2 = GrassmannMPS(scalartype(x), length(lattice2))
	if RealTimeOrderingStyle(lattice2) isa TimeAscending
		# the new time block is appended at the right end, old blocks keep
		# their site positions
		x2.data[1:L] = x.data
	else
		# time descending: the new time block is inserted next to the
		# left boundary, old blocks shift to the right
		x2.data[end-L+1:end] = x.data
	end
	return x2
end

"""
	timesteps(gmps::GrassmannMPS, x::AbstractRealGrassmannLattice)

Return the current time step k of GrassmannMPS
"""
function timesteps(gmps::GrassmannMPS, x::AbstractRealGrassmannLattice)
	L = length(gmps)
	@assert L % (4*x.bands) == 0
	return div(L, 4*x.bands)
end


# a\bar{a}b\bar{b} a₂^+b₂^+ā₂^-b̄₂^-ā₂^+b̄₂^+a₂^-b₂^-  a₁^+b₁^+ā₁^-b̄₁^-ā₁^+b̄₁^+a₁^-b₁^-
function _index(x::AbstractRealGrassmannLattice{<:A1B1ā1b̄1Ā1B̄1a1b1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if branch == :+
		ifelse(conj, TL-4i*bands+2*bands+band, TL-4i*bands+band)
	else
		ifelse(conj, TL-4i*bands+bands+band, TL-4i*bands+3*bands+band)
	end
end

# a\bar{a}b\bar{b} a_2^+\bar{a}_2^+a_2^-\bar{a}_2^-b_2^+\bar{b}_2^+b_2^-\bar{b}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-b_1^+\bar{b}_1^+b_1^-\bar{b}_1^-
function _index(x::AbstractRealGrassmannLattice{<:A1Ā1a1ā1B1B̄1b1b̄1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if branch == :+
		ifelse(conj, TL-4i*bands+2+4*(band-1), TL-4i*bands+1+4*(band-1))
	else
		ifelse(conj, TL-4i*bands+4+4*(band-1), TL-4i*bands+3+4*(band-1))
	end
end

# a\bar{a}b\bar{b} a_1^-\bar{a}_1^-a_1^+\bar{a}_1^+b_1^-\bar{b}_1^-b_1^+\bar{b}_1^+ a_2^-\bar{a}_2^-a_2^+\bar{a}_2^+b_2^-\bar{b}_2^-b_2^+\bar{b}_2^+
# time ascending: the block of time step 1 directly follows the boundary block
function _index(x::AbstractRealGrassmannLattice{<:a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	bands = x.bands
	base = 4*(i-1)*bands
	if branch == :-
		ifelse(conj, base+2+4*(band-1), base+1+4*(band-1))
	else
		ifelse(conj, base+4+4*(band-1), base+3+4*(band-1))
	end
end

# aābb̄ a₂^+ā₂^+b₂^+b̄₂^+a₂^-ā₂^-b₂^-b̄₂^- a₁^+ā₁^+b₁^+b̄₁^+a₁^-ā₁^-b₁^-b̄₁^-
function _index(x::AbstractRealGrassmannLattice{<:A1Ā1B1B̄1a1ā1b1b̄1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if branch == :+
		ifelse(conj, TL-4i*bands+2+2*(band-1), TL-4i*bands+1+2*(band-1))
	else
		ifelse(conj, TL-4i*bands+2*bands+2+2*(band-1), TL-4i*bands+2*bands+1+2*(band-1))
	end
end


function _index(x::AbstractRealGrassmannLattice{<:Ā2A1ā1a2B̄2B1b̄1b̄2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if (i == 1) && conj && (branch==:+)
		length(x) - 2bands + 2*(band-1)+1
	elseif (i == 1) && (!conj) && (branch==:-)
		length(x) - 2bands + 2*band
	elseif (i == x.k) && (!conj) && (branch==:+)
		2*(band-1)+1
	elseif (i == x.k) && conj && (branch==:-)
		2*band
	else
		if branch == :+
			ifelse(conj, (x.k-i)*4*bands + 4*(band-1)+1, (x.k-i-1)*4*bands + 4*(band-1)+2) + 2bands
		else
			ifelse(conj, (x.k-i-1)*4*bands + 4*(band-1)+3, (x.k-i)*4*bands + 4*band) + 2bands
		end
	end
end

# key is timestep, conj, branch, band
function indexmappings(lattice::AbstractRealGrassmannLattice)
	r = Dict{Tuple{Int, Bool, Symbol, Int}, Int}()
	for i in 1:lattice.k
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

function band_boundary(lattice::AbstractRealGrassmannLattice{<:A1Ā1a1ā1B1B̄1b1b̄1}, j::Int)
    posa = lattice[GrassmannContourIndex(j, :+, 1, false)]
    posb = lattice[GrassmannContourIndex(j, :-, lattice.bands, true)]
    return posa, posb
end
function band_boundary(lattice::AbstractRealGrassmannLattice{<:A1Ā1B1B̄1a1ā1b1b̄1}, j::Int)
    posa = lattice[GrassmannContourIndex(j, :+, 1, false)]
    posb = lattice[GrassmannContourIndex(j, :-, lattice.bands, true)]
    return posa, posb
end
function band_boundary(lattice::AbstractRealGrassmannLattice{<:a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2}, j::Int)
    posa = lattice[GrassmannContourIndex(j, :-, 1, false)]
    posb = lattice[GrassmannContourIndex(j, :+, lattice.bands, true)]
    return posa, posb
end
