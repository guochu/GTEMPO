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

"""
	cached_contour_ordered_gf(lattice, a::ContourIndex, b::ContourIndex, A, B...; cache=environments(lattice, A, B...))

Contour-ordered Green's function `G(a, b)` with cached environments. Requires
`conj(a)=false` and `conj(b)=true`; the lesser/greater case is selected by the
ordering of `a` and `b` on the contour.
"""
function cached_contour_ordered_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                                    cache::AbstractExpectationCache=environments(lattice, A, B...)) 
    ((!a.conj) && (b.conj)) || throw(ArgumentError("conj(a)=false and conj(b)=true should be satisfied"))
    return (a < b) ? -cached_gf(lattice, (b, a), A, B...; cache=cache) : cached_gf(lattice, (a, b), A, B...; cache=cache) 
end

"""
	cached_greater(lattice, i::Int, j::Int, A, B...; band=1, kwargs...)

Greater Green's function `G^>(i, j) = -i⟨d(i) d†(j)⟩` on the real branches,
evaluated with cached environments. `band` may be an `Int` or a 2-tuple.
"""
function cached_greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                        band::Union{Int, Tuple{Int, Int}}=1, kwargs...)
    @assert i >= j
    band isa Int && (band = (band, band))
    a = ContourIndex(i, conj=false, branch=:+, band=band[1])
    b = ContourIndex(j, conj=true, branch=:+, band=band[2])
    return cached_gf(lattice, (a, b), A, B...; kwargs...)
end
cached_greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) = cached_greater(lattice, i, 1, A, B...; kwargs...)
"""
	cached_lesser(lattice, i::Int, j::Int, A, B...; band=1, kwargs...)

Lesser Green's function `G^<(i, j) = i⟨d†(i) d(j)⟩` on the real branches,
evaluated with cached environments. `band` may be an `Int` or a 2-tuple.
"""
function cached_lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                       band::Union{Int, Tuple{Int, Int}}=1, kwargs...)
    @assert i <= j
    band isa Int && (band = (band, band))
    a = ContourIndex(i, conj=true, branch=:-, band=band[1])
    b = ContourIndex(j, conj=false, branch=:+, band=band[2])
    return cached_gf(lattice, (a, b), A, B...; kwargs...)
end
cached_lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) = cached_lesser(lattice, 1, i, A, B...; kwargs...)
