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
