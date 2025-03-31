abstract type ImagFockLattice{O<:ImagFockOrdering} <: AbstractFockLattice{O} end
TK.scalartype(::Type{<:ImagFockLattice}) = Float64
branches(::Type{<:ImagFockLattice}) = (:τ,)
TimeOrderingStyle(x::ImagFockLattice) = ImaginaryTimeOrderingStyle(x)

struct ImagFockLattice1Order{O<:ImagFockOrdering} <: ImagFockLattice{O}
	δτ::Float64
	bands::Int
	N::Int
	ordering::O

	ImagFockLattice1Order(δτ::Real, bands::Int, N::Int, ordering::ImagFockOrdering) = new{typeof(ordering)}(convert(Float64, δτ), bands, N, ordering)
end

ImagFockLattice1Order(; δτ::Real, N::Int, bands::Int=1, ordering::ImagFockOrdering=M1N1()) = ImagFockLattice1Order(δτ, bands, N, ordering)
Base.similar(x::ImagFockLattice1Order; δτ::Real=x.δτ, bands::Int=x.bands, N::Int=x.N, ordering::ImagFockOrdering=x.ordering) = ImagFockLattice1Order(δτ, bands, N, ordering)
similargrassmannlattice(x::ImagFockLattice1Order, δτ::Real=x.δτ, bands::Int=x.bands, N::Int=x.N, 
						ordering::ImagGrassmannOrdering=similargrassmannordering(x.ordering)) = GrassmannLattice(contour=:imag, δτ=δτ, N=N, bands=bands, ordering=ordering)

function ImagFockLattice(; order::Int=1, kwargs...)
	(order in (1, 2)) || throw(ArgumentError("order must be 1 or 2"))
	if order == 1
		return ImagFockLattice1Order(; kwargs...)
	else
		error("Second orderr ImagGrassmannLattice not implemented")
	end
end

Base.length(x::ImagFockLattice) = x.bands * x.N
function Base.getproperty(x::ImagFockLattice1Order, s::Symbol)
	if s == :τs
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


function index(x::ImagFockLattice{<:M1N1}, i::Int; band::Int=1, branch::Symbol=:τ)
	@boundscheck begin
		(1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
		(1 <= i <= x.N) || throw(BoundsError(1:x.N, i))
		(branch == :τ) || throw(ArgumentError("branch must be :τ"))
	end
	# TL = length(x)
	return (x.N-i)*x.bands + band
end



function indexmappings(lattice::ImagFockLattice)
	r = Dict{Tuple{Int, Symbol, Int}, Int}()
	for i in 1:lattice.N
		for band in 1:lattice.bands
			r[(i, :τ, band)] = index(lattice, i, band=band)
		end
	end
	return r
end
