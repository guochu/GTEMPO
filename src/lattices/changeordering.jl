"""
    changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, zs::GrassmannMPS...; trunc::TruncationScheme) where {A}

Change the ordering of GVs of the GMPSs x, y, zs from the ordering in lattice to the new ordering A
Inside this function, truncation may not done since swap gate may be used, and "trunc" is used for the
truncation
"""
function changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, zs::Vararg{GrassmannMPS}; 
                        trunc::TruncationScheme=Defaults.alg_trunc_sys()) where {A <: GrassmannOrdering} 
    if OrderingStyle(lattice) isa A
        return lattice, x, y, zs...
    else
        lattice2 = similar(lattice, ordering=A())
        perm = matchindices2(lattice2, lattice)
        return lattice2, permute(x, perm, trunc=trunc), permute(y, perm, trunc=trunc), map(z->permute(z, perm, trunc=trunc), zs)...
    end
end
function changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS; trunc::TruncationScheme=Defaults.alg_trunc_sys()) where {A <: GrassmannOrdering} 
    if OrderingStyle(lattice) isa A
        return lattice, x
    else
        lattice2 = similar(lattice, ordering=A())
        return lattice2, permute(x, matchindices2(lattice2, lattice), trunc=trunc)
    end
end

changeordering(o::GrassmannOrdering, lattice::AbstractGrassmannLattice, x::Vararg{GrassmannMPS}; kwargs...) = changeordering(typeof(o), lattice, x...; kwargs...)

"""
    toadjacentordering(lattice::AbstractGrassmannLattice, x::GrassmannMPS...; kwargs...)

Convert the lattice into a "closest" adjacent ordering
This function a simple wrapper of changeordering, by specifying a particular ordering
"""
toadjacentordering(lattice::AbstractGrassmannLattice, x::Vararg{GrassmannMPS}; kwargs...) = error("toadjacentordering not implemented for lattice type $(typeof(lattice))")
toadjacentordering(lattice::ImagGrassmannLattice, x::Vararg{GrassmannMPS}; kwargs...) = changeordering(A1Ā1B1B̄1, lattice, x...; kwargs...)
toadjacentordering(lattice::RealGrassmannLattice, x::Vararg{GrassmannMPS}; kwargs...) = changeordering(A1Ā1a1ā1B1B̄1b1b̄1, lattice, x...; kwargs...)

