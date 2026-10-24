# local helper: the ⟨aᵢ bⱼ⟩ Green's function built on cached_gf with contour indices
function _cached_gf_ij(lattice::AbstractGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                       b1::Symbol, b2::Symbol, c1::Bool, c2::Bool, band::Union{Int, Tuple{Int, Int}}=1, kwargs...)
    band isa Int && (band = (band, band))
    a = GrassmannContourIndex(i, conj=c1, branch=b1, band=band[1])
    b = GrassmannContourIndex(j, conj=c2, branch=b2, band=band[2])
    return cached_gf(lattice, (a, b), A, B...; kwargs...)
end

"""
	cached_electriccurrent(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, k::Int, A, B...;
                    	               cache=environments(lattice, A, B...), band=1)

Electric current at time step `k`, evaluated as a sum of two-point correlators
with cached environments. Without the `k` argument the whole sequence
`k = 2:lattice.k` is returned.
"""
function cached_electriccurrent(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                	               cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1)
    curr = complex(0.)
    η⁺⁺ = branch(corr; b1=:+, b2=:+)
    η⁺⁻ = branch(corr; b1=:+, b2=:-)
    for j in 1:k-1
        curr += η⁺⁺[k, j] * _cached_gf_ij(lattice, k, j, A, B...; cache=cache, b1=:+, b2=:+, band=band, c1=true, c2=false)
        curr += η⁺⁻[k, j] * _cached_gf_ij(lattice, k, j, A, B...; cache=cache, b1=:+, b2=:-, band=band, c1=true, c2=false)
    end
    return 2 * curr / lattice.δt
end
cached_electriccurrent(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                        cache::AbstractExpectationCache=environments(lattice, A, B...), 
                        kwargs...) = [cached_electriccurrent(lattice, corr, k, A, B...; cache=cache, kwargs...) for k in 2:lattice.k]

"""
	cached_electriccurrent_fast(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, k::Int, A, B...;
                                    cache=environments(lattice, A, B...), band=1)

Electric current at time step `k` computed with the MPO-based (fast) current
evaluation. Without the `k` argument the whole sequence is returned.
"""
function cached_electriccurrent_fast(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                                    cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1)
    mpo = build_current_mpo(lattice, corr, k, band)
    curr = expectationvalue(mpo, cache)
    return 2 * curr / lattice.δt
end
cached_electriccurrent_fast(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                        cache::AbstractExpectationCache=environments(lattice, A, B...), 
                        kwargs...) = [cached_electriccurrent_fast(lattice, corr, k, A, B...; cache=cache, kwargs...) for k in 2:lattice.k]


"""
	cached_heatcurrent_fast(lattice::RealGrassmannLattice, bath::AbstractFermionicNormalBath, args...; kwargs...)

Heat current computed with the cached MPO-based (fast) current evaluation,
using the `ω`-weighted correlation function built by `heatcorrelationfunction`.
"""
cached_heatcurrent_fast(lattice::RealGrassmannLattice, bath::AbstractFermionicNormalBath, args...; kwargs...) = cached_electriccurrent_fast(
                        lattice, heatcorrelationfunction(bath, lattice), args...; kwargs...)



