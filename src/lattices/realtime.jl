abstract type RealGrassmannLattice{O<:RealGrassmannOrdering} <: AbstractGrassmannLattice{O} end
TK.scalartype(::Type{<:RealGrassmannLattice}) = ComplexF64
branches(::Type{<:RealGrassmannLattice}) = (:+, :-)

# k is the number of discretization, nbands is the number of bands
# pos is the position within a band
# k+1 due to the grassmann number on the boundary for the final trace

"""
	struct RealGrassmannLattice1Order <: RealGrassmannLattice

First order splitting of the real-time contour
"""
struct RealGrassmannLattice1Order{O<:RealGrassmannOrdering} <: RealGrassmannLattice{O}
	δt::Float64
	bands::Int
	N::Int
	ordering::O

	RealGrassmannLattice1Order(δt::Real, bands::Int, N::Int, ordering::RealGrassmannOrdering) = new{typeof(ordering)}(convert(Float64, δt), bands, N, ordering)
end

# the default is that the system starts from 0 temperature (state 0)
RealGrassmannLattice1Order(; δt::Real, N::Int, bands::Int=1, ordering::RealGrassmannOrdering=A1A1a1a1B1B1b1b1()) = RealGrassmannLattice1Order(δt, bands, N, ordering)
Base.similar(x::RealGrassmannLattice1Order; δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N, ordering::RealGrassmannOrdering=x.ordering) = RealGrassmannLattice1Order(δt, bands, N, ordering)

"""
	struct RealGrassmannLattice2Order <: RealGrassmannLattice

Second order splitting of the real-time contour,
we use a symmetric splitting with the same number
of discretization as first order, however, the 
MPS-IF has to grow during runtime
"""
struct RealGrassmannLattice2Order{O<:RealGrassmannOrdering} <: RealGrassmannLattice{O}
	δt::Float64
	bands::Int
	N::Int
	ordering::O

	RealGrassmannLattice2Order(δt::Real, bands::Int, N::Int, ordering::RealGrassmannOrdering) = new{typeof(ordering)}(convert(Float64, δt), bands, N, ordering)
end
RealGrassmannLattice2Order(; δt::Real, N::Int=0, bands::Int=1, ordering::RealGrassmannOrdering=A1A1a1a1B1B1b1b1()) = RealGrassmannLattice2Order(δt, bands, N, ordering)
Base.similar(x::RealGrassmannLattice2Order; δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N, ordering::RealGrassmannOrdering=x.ordering) = RealGrassmannLattice2Order(δt, bands, N, ordering)

function Base.getproperty(x::RealGrassmannLattice, s::Symbol)
	if (s == :k) || (s == :kt)
		return x.N+1
	elseif s == :t
		return x.N * x.δt
	elseif s == :ts
		return 0:x.δt:x.N*x.δt
	elseif s == :Nt
		return x.N
	else
		getfield(x, s)
	end
end
Base.length(x::RealGrassmannLattice) = 4*x.bands * x.k + 2 * x.bands
# Base.length(x::RealGrassmannLattice) = _length(x.bands, x.k)
# _length(bands::Int, k::Int) = 4*bands*k + 2*bands
function RealGrassmannLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return RealGrassmannLattice1Order(; kwargs...)
	else
		# error("Second orderr RealGrassmannLattice not implemented")
		return RealGrassmannLattice2Order(; kwargs...)
	end
end

makestep(x::RealGrassmannLattice1Order) = RealGrassmannLattice1Order(x.δt, x.bands, x.N+1, x.ordering)
makestep(x::RealGrassmannLattice2Order) = RealGrassmannLattice2Order(x.δt, x.bands, x.N+1, x.ordering)


"""
	makestep(x::GrassmannMPS, lattice::RealGrassmannLattice)
Increase the time step of MPS-IF and the lattice by 1
"""
function makestep(lattice::RealGrassmannLattice, x::GrassmannMPS...)
	(all(v->length(v)==length(lattice), x)) || throw(DimensionMismatch())
	(LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("makestep only implemented for TimeLocalLayout"))
	# (length(x) < length(lattice)) || throw(ArgumentError("MPS-IF has reached amximum length"))
	lattice2 = makestep(lattice)
	x2 = map(y->_makestep_util(lattice2, y), x)
	return lattice2, x2...
end

function _makestep_util(lattice2::RealGrassmannLattice, x::GrassmannMPS)
	L = length(x)
	x2 = GrassmannMPS(scalartype(x), length(lattice2))
	x2.data[end-L+1:end] = x.data
	return x2
end

"""
	timesteps(gmps::GrassmannMPS, x::RealGrassmannLattice)

Return the current time step k of GrassmannMPS
"""
function timesteps(gmps::GrassmannMPS, x::RealGrassmannLattice)
	L = length(gmps) - 2 * x.bands
	@assert L % (4*x.bands) == 0
	return div(L, 4*x.bands)
end

# ab\bar{b}\bar{a} a_2^+a_2^-b_2^+b_2^-\bar{b}_2^-\bar{b}_2^+\bar{a}_2^-\bar{a}_2^+ a_1^+a_1^-b_1^+b_1^-\bar{b}_1^-\bar{b}_1^+\bar{a}_1^-\bar{a}_1^+
function index(x::RealGrassmannLattice{<:A1a1B1b1b1B1a1A1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2bands+1-band, band)
	else
		if branch == :+
			ifelse(conj, TL-4(i-1)*bands-2(band-1), TL-4i*bands+1+2(band-1))
		else
			ifelse(conj, TL-4(i-1)*bands-1-2(band-1), TL-4i*bands+2band)
		end
	end
end

# a\bar{a}b\bar{b} a₂^+b₂^+ā₂^-b̄₂^-ā₂^+b̄₂^+a₂^-b₂^-  a₁^+b₁^+ā₁^-b̄₁^-ā₁^+b̄₁^+a₁^-b₁^-
function index(x::RealGrassmannLattice{<:A1B1ā1b̄1A1B1a1b1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1)
	else
		if branch == :+
			ifelse(conj, TL-4i*bands+2*bands+band, TL-4i*bands+band)
		else
			ifelse(conj, TL-4i*bands+bands+band, TL-4i*bands+3*bands+band)
		end
	end
end

# a\bar{a}b\bar{b} a_2^+\bar{a}_2^+a_2^-\bar{a}_2^-b_2^+\bar{b}_2^+b_2^-\bar{b}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-b_1^+\bar{b}_1^+b_1^-\bar{b}_1^-
function index(x::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1)
	else
		if branch == :+
			ifelse(conj, TL-4i*bands+2+4*(band-1), TL-4i*bands+1+4*(band-1))
		else
			ifelse(conj, TL-4i*bands+4+4*(band-1), TL-4i*bands+3+4*(band-1))
		end
	end
end

# aābb̄ a₂^+ā₂^+b₂^+b̄₂^+a₂^-ā₂^-b₂^-b̄₂^- a₁^+ā₁^+b₁^+b̄₁^+a₁^-ā₁^-b₁^-b̄₁^-
function index(x::RealGrassmannLattice{<:A1A1B1B1a1a1b1b1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1)
	else
		if branch == :+
			ifelse(conj, TL-4i*bands+2+2*(band-1), TL-4i*bands+1+2*(band-1))
		else
			ifelse(conj, TL-4i*bands+2*bands+2+2*(band-1), TL-4i*bands+2*bands+1+2*(band-1))
		end
	end
end

# a\bar{a}b\bar{b} a_2^+\bar{a}_2^+a_1^+\bar{a}_1^+ a_2^-\bar{a}_2^-a_1^-\bar{a}_1^- b_2^+\bar{b}_2^+b_1^+\bar{b}_1^+  b_2^-\bar{b}_2^-b_1^-\bar{b}_1^-
function index(x::RealGrassmannLattice{<:A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	n = 4 * x.k
	if i == 0
		ifelse(conj, 2*band, 2*band-1)
	else
		if branch == :+
			ifelse(conj, (band-1) * n + 2*(x.k-i) + 2, (band-1) * n + 2*(x.k-i) + 1 ) + 2*x.bands
		else
			ifelse(conj, (band-1) * n + 2*x.k+ 2*(x.k-i) + 2, (band-1) * n + 2*x.k+ 2*(x.k-i) + 1 ) + 2*x.bands
		end
	end
end

# ab\bar{b}\bar{a} a_2^+b_2^+\bar{b}_2^+\bar{a}_2^+a_1^+b_1^+\bar{b}_1^+\bar{a}_1^+ a_1^-b_1^-\bar{b}_1^-\bar{a}_1^-a_2^-b_2^-\bar{b}_2^-\bar{a}_2^-
function index(x::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2*bands+1-band, band)
	else
		if branch == :+
			ifelse(conj, 2*bands*(x.k-i) + 2*bands-band+1, 2*bands*(x.k-i) + band ) + 2*bands
		else
			ifelse(conj, 2*bands*(i-1) + 2*bands-band+1, 2*bands*(i-1) + band ) + 2*bands + 2 * bands * x.k
		end
	end
end


# a\bar{a}b\bar{b} a_2^+\bar{a}_2^+b_2^+\bar{b}_2^+a_1^+\bar{a}_1^+b_1^+\bar{b}_1^+ a_1^-\bar{a}_1^-b_1^-\bar{b}_1^-a_2^-\bar{a}_2^-b_2^-\bar{b}_2^-
function index(x::RealGrassmannLattice{<:A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(0 <= i <= x.k) || throw(BoundsError(0:x.k, i))
		(branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
	end
	TL = length(x)
	bands = x.bands
	if i == 0
		ifelse(conj, 2*band, 2*band-1)
	else
		if branch == :+
			ifelse(conj, 2*bands*(x.k-i) + 2*band, 2*bands*(x.k-i) + 2*band - 1 ) + 2*bands
		else
			ifelse(conj, 2*bands*(i-1) + 2*band, 2*bands*(i-1) + 2*band - 1 ) + 2*bands + 2 * bands * x.k
		end
	end
end

# key is timestep, conj, branch, band
function indexmappings(lattice::RealGrassmannLattice)
	r = Dict{Tuple{Int, Bool, Symbol, Int}, Int}()
	for i in 0:lattice.k
		for c in (true, false)
			for band in 1:lattice.bands
				if i == 0
					r[(i, c, :+, band)] = index(lattice, i, conj=c, band=band)
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

function band_boundary(lattice::RealGrassmannLattice{<:A1a1B1b1b1B1a1A1}, j::Int)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=1)     
    else
        posa = index(lattice, j, conj=false, branch=:+, band=1)
        posb = index(lattice, j, conj=true, branch=:+, band=1)
    end
    return posa, posb
end
function band_boundary(lattice::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, j::Int)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=lattice.bands)     
    else
        posa = index(lattice, j, conj=false, branch=:+, band=1)
        posb = index(lattice, j, conj=true, branch=:-, band=lattice.bands)
    end
    return posa, posb
end
function band_boundary(lattice::RealGrassmannLattice{<:A1A1B1B1a1a1b1b1}, j::Int)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=lattice.bands)     
    else
        posa = index(lattice, j, conj=false, branch=:+, band=1)
        posb = index(lattice, j, conj=true, branch=:-, band=lattice.bands)
    end
    return posa, posb
end
function band_boundary(lattice::RealGrassmannLattice{<:A1B1ā1b̄1A1B1a1b1}, j::Int)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=lattice.bands)     
    else
        posa = index(lattice, j, conj=false, branch=:+, band=1)
        posb = index(lattice, j, conj=false, branch=:-, band=lattice.bands)
    end
    return posa, posb
end
function band_boundary(lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, j::Int; branch::Symbol=:+)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=1)   
    else
        posa = index(lattice, j, conj=false, band=1, branch=branch)
        posb = index(lattice, j, conj=true, band=1, branch=branch) 
    end
    return posa, posb
end
function band_boundary(lattice::RealGrassmannLattice{<:A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, j::Int; branch::Symbol=:+)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=lattice.bands)   
    else
        posa = index(lattice, j, conj=false, band=1, branch=branch)
        posb = index(lattice, j, conj=true, band=lattice.bands, branch=branch) 
    end
    return posa, posb
end
