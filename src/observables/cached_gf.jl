# cached version of calculating Green's functions

"""
   cached_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...)

The same as gf, but use a cache for efficiency, the cache can be precomputed
and used for calculating any other observables
"""
function cached_gf(lattice::AbstractGrassmannLattice, a::NTuple{N, ContourIndex}, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...)) where {N}
    pos = map(x->lattice[x], a)
    t = GTerm(pos, coeff=1)
    return expectationvalue(t, cache)
end


# function cached_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
#                     cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...)
#     pos1, pos2 = lattice[a], lattice[b]
#     t = GTerm(pos1, pos2, coeff=1)
#     return expectationvalue(t, cache; kwargs...)
# end

function cached_contour_ordered_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                                    cache::AbstractExpectationCache=environments(lattice, A, B...)) 
    ((!a.conj) && (b.conj)) || throw(ArgumentError("conj(a)=false and conj(b)=true should be satisfied"))
    return (a < b) ? -cached_gf(lattice, (b, a), A, B...; cache=cache) : cached_gf(lattice, (a, b), A, B...; cache=cache) 
end

# imaginary time
function cached_Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), c1::Bool=false, c2::Bool=true, band::Union{Int, Tuple{Int, Int}}=1)
    if isa(band, Int)
        band = (band, band)
    end
    band1, band2 = band    
	a, b = ContourIndex(i, conj=c1, branch=:τ, band=band1), ContourIndex(j, conj=c2, branch=:τ, band=band2)
    return cached_gf(lattice, (a, b), A, B...; cache=cache)
end
cached_Gτ(lattice::ImagGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) = cached_Gτ(lattice, i, 1, A, B...; kwargs...)


# real time
function cached_Gt(lattice::RealGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol, b2::Symbol, c1::Bool=true, 
                    c2::Bool=false, band::Union{Int, Tuple{Int, Int}}=1)
    if isa(band, Int)
        band = (band, band)
    end
    band1, band2 = band    
    a, b = ContourIndex(i, conj=c1, branch=b1, band=band1), ContourIndex(j, conj=c2, branch=b2, band=band2)
    return cached_gf(lattice, (a, b), A, B...; cache=cache)
end

function cached_Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol, b2::Symbol, c1::Bool=true, 
                    c2::Bool=false, band::Union{Int, Tuple{Int, Int}}=1)
    if isa(band, Int)
        band = (band, band)
    end
    band1, band2 = band    
    a, b = ContourIndex(i, conj=c1, branch=b1, band=band1), ContourIndex(j, conj=c2, branch=b2, band=band2)
    return cached_gf(lattice, (a, b), A, B...; cache=cache)
end

function cached_contour_ordered_Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol, b2::Symbol, band::Union{Int, Tuple{Int, Int}}=1) 
    if isa(band, Int)
        band = (band, band)
    end
    band1, band2 = band    
    a, b = ContourIndex(i, conj=false, branch=b1, band=band1), ContourIndex(j, conj=true, branch=b2, band=band2)
    return cached_contour_ordered_gf(lattice, a, b, A, B...; cache=cache)
end

#########*******************************###########

function cached_Gτ(lattice::Union{ImagGrassmannLattice, MixedGrassmannLattice}, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                    cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, kwargs...)
	g = zeros(_scalartype(A), lattice.kτ)
	for i in 1:lattice.kτ-1
		g[i] = cached_Gτ(lattice, i, A, B...; cache=cache, band=band, kwargs...)
	end
	g[end] = 1 - g[1]
	return g	
end


function cached_Gt(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; b1::Symbol, b2::Symbol, kwargs...)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    return cached_Gm(lattice, i, j, A, B...; b1=b1, b2=b2, kwargs...)
end
function cached_Gτ(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; c1::Bool=false, c2::Bool=true, kwargs...)
    return cached_Gm(lattice, i, j, A, B...; b1=:τ, b2=:τ, c1=c1, c2=c2, kwargs...)
end
cached_Gτ(lattice::MixedGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) = cached_Gτ(lattice, i, 1, A, B...; kwargs...)

function cached_greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...)
    @assert i >= j
    return cached_Gt(lattice, i, j, A, B...; b1=:+, b2=:+, c1=false, c2=true, kwargs...)
end
cached_greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) = cached_greater(lattice, i, 1, A, B...; kwargs...)
function cached_lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...)
    @assert i <= j
    return cached_Gt(lattice, i, j, A, B...; b1=:-, b2=:+, c1=true, c2=false, kwargs...)
end
cached_lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) = cached_lesser(lattice, 1, i, A, B...; kwargs...)

