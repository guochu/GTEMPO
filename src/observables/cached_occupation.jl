# real-time first order
function cached_occupation(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) 
    return real(cached_Gt(lattice, i, i, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end
cached_occupation(lattice::RealGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) = [
                    cached_occupation(lattice, i, A, B...; cache=cache, kwargs...) for i in 1:lattice.N]
cached_occupation(lattice::MixedGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) = [
                    cached_occupation(lattice, i, A, B...; cache=cache, kwargs...) for i in 1:lattice.Nt]


# real-time second order
function cached_occupation(lattice::RealGrassmannLattice2Order, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; kwargs...) 
    return real(cached_Gt(lattice, lattice.k, lattice.k, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end