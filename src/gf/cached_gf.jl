
# imaginary time
function cached_gf(lattice::ImagGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, kwargs...)
	pos1, pos2 = index(lattice, i, conj=false, band=band), index(lattice, 1, conj=true, band=band)
	t = GTerm(pos1, pos2, coeff=1)
	return expectationvalue(t, cache; kwargs...)
end

function cached_gf(lattice::ImagGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...;
                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...)
	g = zeros(Float64, lattice.k)
	for i in 1:lattice.k-1
		g[i] = cached_gf(lattice, i, A, B...; cache=cache, kwargs...)
	end
	g[end] = 1 - g[1]
	return g	
end

# real time
function cached_gf(lattice::RealGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Int=1, kwargs...)
    pos1, pos2 = index(lattice, i, conj=c1, branch=b1, band=band), index(lattice, j, conj=c2, branch=b2, band=band)
    t = GTerm(pos1, pos2, coeff=1)
    return expectationvalue(t, cache; kwargs...)
end

# real-time first order
function cached_occupation(lattice::RealGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) 
    return real(cached_gf(lattice, i, i, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end
cached_occupation(lattice::RealGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) = [
        cached_occupation(lattice, i, A, B...; cache=cache, kwargs...) for i in 1:lattice.N]


function cached_electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                	               cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, max_range::Int=10000, kwargs...)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    for j in max(1, k-max_range):k-1
        curr += η⁺⁺[k, j] * cached_gf(lattice, k, j, A, B...; cache=cache, b1=:+, b2=:+, band=band, kwargs...)
        curr += η⁺⁻[k, j] * cached_gf(lattice, k, j, A, B...; cache=cache, b1=:+, b2=:-, band=band, kwargs...)
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


# real-time second order
function cached_occupation(lattice::RealGrassmannLattice2Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) 
    return real(cached_gf(lattice, lattice.k, lattice.k, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end
function cached_electriccurrent(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                                cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, max_range::Int=10000, kwargs...)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    @assert  lattice.N <= div(size(η⁺⁺,1)-1, 2) 
    k = lattice.k
    curr -= η⁺⁺[2*k-1, 1] * cached_gf(lattice, k, 1, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
    curr -= η⁺⁻[2*k-1, 1] * cached_gf(lattice, k, 1, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)
    for j in 2:k-1
        curr -= (η⁺⁺[2*k-1, 2*j-2] + η⁺⁺[2*k-1, 2*j-1]) * cached_gf(lattice, k, j, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
        curr -= (η⁺⁻[2*k-1, 2*j-2] + η⁺⁻[2*k-1, 2*j-1])  * cached_gf(lattice, k, j, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)
    end
    curr -= η⁺⁺[2*k-1, 2*k-2] * cached_gf(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
    curr -= η⁺⁻[2*k-1, 2*k-2] * cached_gf(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)
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
    curr -= η⁺⁺[2*k-1, 2*k-2] * cached_gf(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
    curr -= η⁺⁻[2*k-1, 2*k-2] * cached_gf(lattice, k, k, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)

    return 2 * curr / (0.5*lattice.δt)
end
