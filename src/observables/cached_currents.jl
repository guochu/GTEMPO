function cached_electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                	               cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    for j in 1:k-1
        curr += η⁺⁺[k, j] * cached_Gt(lattice, k, j, A, B...; cache=cache, b1=:+, b2=:+, band=band, c1=true, c2=false)
        curr += η⁺⁻[k, j] * cached_Gt(lattice, k, j, A, B...; cache=cache, b1=:+, b2=:-, band=band, c1=true, c2=false)
    end
    return 2 * curr / lattice.δt
end
cached_electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                        cache::AbstractExpectationCache=environments(lattice, A, B...), 
                        kwargs...) = [cached_electriccurrent(lattice, corr, k, A, B...; cache=cache, kwargs...) for k in 2:lattice.k]

function cached_electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                                    cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1)
    mpo = build_current_mpo(lattice, corr, k, band)
    curr = expectationvalue(mpo, cache)
    return 2 * curr / lattice.δt
end
cached_electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                        cache::AbstractExpectationCache=environments(lattice, A, B...), 
                        kwargs...) = [cached_electriccurrent_fast(lattice, corr, k, A, B...; cache=cache, kwargs...) for k in 2:lattice.k]


cached_heatcurrent_fast(lattice::RealGrassmannLattice1Order, bath::AbstractFermionicNormalBath, args...; kwargs...) = cached_electriccurrent_fast(
                        lattice, heatcorrelationfunction(bath, lattice), args...; kwargs...)


function cached_electriccurrent(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                                cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    @assert  lattice.N <= div(size(η⁺⁺,1)-1, 2) 
    k = lattice.k
    curr -= η⁺⁺[2*k-1, 1] * cached_Gt(lattice, k, 1, A, B...; cache=cache, b1=:-, b2=:+, band=band, c1=true, c2=false)
    curr -= η⁺⁻[2*k-1, 1] * cached_Gt(lattice, k, 1, A, B...; cache=cache, b1=:-, b2=:-, band=band, c1=true, c2=false)
    for j in 2:k-1
        curr -= (η⁺⁺[2*k-1, 2*j-2] + η⁺⁺[2*k-1, 2*j-1]) * cached_Gt(lattice, k, j, A, B...; cache=cache, b1=:-, b2=:+, band=band, c1=true, c2=false)
        curr -= (η⁺⁻[2*k-1, 2*j-2] + η⁺⁻[2*k-1, 2*j-1])  * cached_Gt(lattice, k, j, A, B...; cache=cache, b1=:-, b2=:-, band=band, c1=true, c2=false)
    end
    curr -= η⁺⁺[2*k-1, 2*k-2] * cached_Gt(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:+, band=band, c1=true, c2=false)
    curr -= η⁺⁻[2*k-1, 2*k-2] * cached_Gt(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:-, band=band, c1=true, c2=false)
    return 2 * curr / (0.5*lattice.δt)
end

function cached_electriccurrent_fast(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                                    cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1)
    k = lattice.k
    mpo = build_current_mpo(lattice, corr, k, band)
    # A2 =_mult_A(mpo, A) 
    # curr = cached_integrate_util(lattice, k, 1, cache, A2, B...; kwargs...)
    # a, b = first(positions(mpo)), last(positions(mpo))
    curr = expectationvalue(mpo, cache)

    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    curr -= η⁺⁺[2*k-1, 2*k-2] * cached_Gt(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:+, band=band, c1=true, c2=false)
    curr -= η⁺⁻[2*k-1, 2*k-2] * cached_Gt(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:-, band=band, c1=true, c2=false)

    return 2 * curr / (0.5*lattice.δt)
end
