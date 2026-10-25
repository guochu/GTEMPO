# real-time first order
"""
	cached_occupation(lattice, i::Int, A, B...; band=1, kwargs...)

Occupation `⟨n̂(i)⟩` at time step `i` on real-time lattices, evaluated with
cached environments. `band` may be an `Int` or a 2-tuple.
"""
function cached_occupation(lattice::Union{AbstractRealGrassmannLattice, AbstractMixedGrassmannLattice}, i::Int, A::GrassmannMPS, B::Vararg{GrassmannMPS}; kwargs...) 
    band = get(kwargs, :band, 1)
    band isa Int && (band = (band, band))
    a = GrassmannContourIndex(i, conj=false, branch=:+, band=band[1])
    b = GrassmannContourIndex(i, conj=true, branch=:-, band=band[2])
    return real(cached_gf(lattice, (a, b), A, B...; kwargs...))
end
cached_occupation(lattice::AbstractRealGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) = [
                    cached_occupation(lattice, i, A, B...; cache=cache, kwargs...) for i in 1:lattice.N]
cached_occupation(lattice::AbstractMixedGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), kwargs...) = [
                    cached_occupation(lattice, i, A, B...; cache=cache, kwargs...) for i in 1:lattice.Nt]
