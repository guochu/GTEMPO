# cached version of calculating Green's functions

"""
   cached_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...)

The same as gf, but use a cache for efficiency, the cache can be precomputed
and used for calculating any other observables
"""
function cached_gf(lattice::AbstractGrassmannLattice, a::NTuple{N, ContourIndex}, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) where {N}
    pos = map(x->lattice[x], a)
    t = GTerm(pos, coeff=1)
    return expectationvalue(t, cache; kwargs...)
end


# function cached_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
#                     cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...)
#     pos1, pos2 = lattice[a], lattice[b]
#     t = GTerm(pos1, pos2, coeff=1)
#     return expectationvalue(t, cache; kwargs...)
# end

function cached_contour_ordered_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) 
    ((!a.conj) && (b.conj)) || throw(ArgumentError("conj(a)=false and conj(b)=true should be satisfied"))
    return (a < b) ? -cached_gf(lattice, (b, a), A, B...; cache=cache, kwargs...) : cached_gf(lattice, (a, b), A, B...; cache=cache, kwargs...) 
end

# imaginary time
function cached_Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                cache::AbstractExpectationCache=environments(lattice, A, B...), c1::Bool=false, c2::Bool=true, band::Int=1, kwargs...)
	a, b = ContourIndex(i, conj=c1, branch=:τ, band=band), ContourIndex(j, conj=c2, branch=:τ, band=band)
    return cached_gf(lattice, (a, b), A, B...; cache=cache, kwargs...)
end
cached_Gτ(lattice::ImagGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) = cached_Gτ(lattice, i, 1, A, B...; kwargs...)


# real time
function cached_Gt(lattice::RealGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Int=1, kwargs...)
    a, b = ContourIndex(i, conj=c1, branch=b1, band=band), ContourIndex(j, conj=c2, branch=b2, band=band)
    return cached_gf(lattice, (a, b), A, B...; cache=cache, kwargs...)
end

function cached_Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Int=1, kwargs...)
    a, b = ContourIndex(i, conj=c1, branch=b1, band=band), ContourIndex(j, conj=c2, branch=b2, band=band)
    return cached_gf(lattice, (a, b), A, B...; cache=cache, kwargs...)
end

function cached_contour_ordered_Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol, b2::Symbol, band::Int=1, kwargs...) 
    a, b = ContourIndex(i, conj=false, branch=b1, band=band), ContourIndex(j, conj=true, branch=b2, band=band)
    return cached_contour_ordered_gf(lattice, a, b, A, B...; cache=cache, kwargs...)
end

#########*******************************###########

function cached_Gτ(lattice::Union{ImagGrassmannLattice, MixedGrassmannLattice}, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...;
                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...)
	g = zeros(scalartype(lattice), lattice.kτ)
	for i in 1:lattice.kτ-1
		g[i] = cached_Gτ(lattice, i, A, B...; cache=cache, kwargs...)
	end
	g[end] = 1 - g[1]
	return g	
end


function cached_Gt(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; b1::Symbol, b2::Symbol, kwargs...)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    return cached_Gm(lattice, i, j, A, B...; b1=b1, b2=b2, kwargs...)
end
function cached_Gτ(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; c1::Bool=false, c2::Bool=true, kwargs...)
    return cached_Gm(lattice, i, j, A, B...; b1=:τ, b2=:τ, c1=c1, c2=c2, kwargs...)
end
cached_Gτ(lattice::MixedGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = cached_Gτ(lattice, i, 1, A, B...; kwargs...)

function cached_greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...)
    @assert i >= j
    return cached_Gt(lattice, i, j, A, B...; b1=:+, b2=:+, c1=false, c2=true, kwargs...)
end
cached_greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = cached_greater(lattice, i, 1, A, B...; kwargs...)
function cached_lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...)
    @assert i <= j
    return cached_Gt(lattice, i, j, A, B...; b1=:-, b2=:+, c1=true, c2=false, kwargs...)
end
cached_lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = cached_lesser(lattice, 1, i, A, B...; kwargs...)


# real-time first order
function cached_occupation(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) 
    return real(cached_Gt(lattice, i, i, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end
cached_occupation(lattice::RealGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) = [
                    cached_occupation(lattice, i, A, B...; cache=cache, kwargs...) for i in 1:lattice.N]
cached_occupation(lattice::MixedGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) = [
                    cached_occupation(lattice, i, A, B...; cache=cache, kwargs...) for i in 1:lattice.Nt]


function cached_electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                	               cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, max_range::Int=10000, kwargs...)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    for j in max(1, k-max_range):k-1
        curr += η⁺⁺[k, j] * cached_Gt(lattice, k, j, A, B...; cache=cache, b1=:+, b2=:+, band=band, kwargs...)
        curr += η⁺⁻[k, j] * cached_Gt(lattice, k, j, A, B...; cache=cache, b1=:+, b2=:-, band=band, kwargs...)
    end
    return 2 * curr / lattice.δt
end
cached_electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                        cache::AbstractExpectationCache=environments(lattice, A, B...), 
                        kwargs...) = [cached_electriccurrent(lattice, corr, k, A, B...; cache=cache, kwargs...) for k in 2:lattice.k]

function cached_electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                                    cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, kwargs...)
    mpo = build_current_mpo(lattice, corr, k, band)
    curr = expectationvalue(mpo, cache; kwargs...)
    return 2 * curr / lattice.δt
end
cached_electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                        cache::AbstractExpectationCache=environments(lattice, A, B...), 
                        kwargs...) = [cached_electriccurrent_fast(lattice, corr, k, A, B...; cache=cache, kwargs...) for k in 2:lattice.k]

function cached_heatcurrent_fast(lattice::RealGrassmannLattice1Order, bath::AbstractFermionicBath, args...; kwargs...)
    bath2 = similar(bath, _mult_w(bath.spectrum))
    corr = correlation(bath2, lattice)
    return cached_electriccurrent_fast(lattice, corr, args...; kwargs...)
end

# real-time second order
function cached_occupation(lattice::RealGrassmannLattice2Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) 
    return real(cached_Gt(lattice, lattice.k, lattice.k, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end
function cached_electriccurrent(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                                cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, max_range::Int=10000, kwargs...)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    @assert  lattice.N <= div(size(η⁺⁺,1)-1, 2) 
    k = lattice.k
    curr -= η⁺⁺[2*k-1, 1] * cached_Gt(lattice, k, 1, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
    curr -= η⁺⁻[2*k-1, 1] * cached_Gt(lattice, k, 1, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)
    for j in 2:k-1
        curr -= (η⁺⁺[2*k-1, 2*j-2] + η⁺⁺[2*k-1, 2*j-1]) * cached_Gt(lattice, k, j, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
        curr -= (η⁺⁻[2*k-1, 2*j-2] + η⁺⁻[2*k-1, 2*j-1])  * cached_Gt(lattice, k, j, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)
    end
    curr -= η⁺⁺[2*k-1, 2*k-2] * cached_Gt(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
    curr -= η⁺⁻[2*k-1, 2*k-2] * cached_Gt(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)
    return 2 * curr / (0.5*lattice.δt)
end

function cached_electriccurrent_fast(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                                    cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, kwargs...)
    k = lattice.k
    mpo = build_current_mpo(lattice, corr, k, band)
    # A2 =_mult_A(mpo, A) 
    # curr = cached_integrate_util(lattice, k, 1, cache, A2, B...; kwargs...)
    # a, b = first(positions(mpo)), last(positions(mpo))
    curr = expectationvalue(mpo, cache; kwargs...)

    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    curr -= η⁺⁺[2*k-1, 2*k-2] * cached_Gt(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
    curr -= η⁺⁻[2*k-1, 2*k-2] * cached_Gt(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)

    return 2 * curr / (0.5*lattice.δt)
end
