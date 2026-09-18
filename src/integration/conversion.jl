# grassmann ordering conversion
# changeordering(::Type{A}, lattice::AbstractGrassmannLattice{A}, x::GrassmannMPS, ys::GrassmannMPS...; kwargs...) where A = x, ys...
# function changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, zs::GrassmannMPS...; kwargs...) where A
#     x1 = convert_ordering(A, lattice, x; kwargs...)
#     y1 = convert_ordering(A, lattice, y; kwargs...)
#     zs = map(z -> convert_ordering(A, lattice, z; kwargs...), zs)
#     return x1, y1, zs...
# end 
# changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS; kwargs...) where A = convert_ordering(A, lattice, x; kwargs...)


"""
    changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, zs::GrassmannMPS...; trunc::TruncationScheme) where {A}

Change the ordering of GVs of the GMPSs x, y, zs from the ordering in lattice to the new ordering A
Inside this function, truncation may not done since swap gate may be used, and "trunc" is used for the
truncation
"""
function changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, zs::Vararg{GrassmannMPS}; 
                        trunc::TruncationScheme=DefaultKTruncation) where {A <: GrassmannOrdering} 
    if OrderingStyle(lattice) isa A
        return lattice, x, y, zs...
    else
        lattice2 = similar(lattice, ordering=A())
        perm = matchindices2(lattice2, lattice)
        return lattice2, permute(x, perm, trunc=trunc), permute(y, perm, trunc=trunc), map(z->permute(z, perm, trunc=trunc), zs)...
    end
end
function changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation) where {A <: GrassmannOrdering} 
    if OrderingStyle(lattice) isa A
        return lattice, x
    else
        lattice2 = similar(lattice, ordering=A())
        return lattice2, permute(x, matchindices2(lattice2, lattice), trunc=trunc)
    end
end

changeordering(o::GrassmannOrdering, lattice::AbstractGrassmannLattice, x::Vararg{GrassmannMPS}; kwargs...) = changeordering(typeof(o), lattice, x...; kwargs...)

# function convert_ordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation) where {A <: GrassmannOrdering} 
#     lattice2 = similar(lattice, ordering=A())
#     return permute(x, matchindices2(lattice2, lattice), trunc=trunc)
# end


# function convert_ordering(::Type{A1Ā1B1B̄1}, lattice::ImagGrassmannLattice{A1B1B̄1Ā1}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
#     return _abba2aabb(x, lattice, trunc=trunc)
# end

#     return _abba2aabb(x, lattice, trunc=trunc)
# end


"""
    toadjacentordering(lattice::AbstractGrassmannLattice, x::GrassmannMPS...; kwargs...)

Convert the lattice into a "closest" adjacent ordering
This function a simple wrapper of changeordering, by specifying a particular ordering 
"""
toadjacentordering(lattice::AbstractGrassmannLattice, x::Vararg{GrassmannMPS}; kwargs...) = error("toadjacentordering not implemented for lattice type $(typeof(lattice))")
toadjacentordering(lattice::ImagGrassmannLattice, x::Vararg{GrassmannMPS}; kwargs...) = changeordering(A1Ā1B1B̄1, lattice, x...; kwargs...)
toadjacentordering(lattice::RealGrassmannLattice, x::Vararg{GrassmannMPS}; kwargs...) = changeordering(A1Ā1a1ā1B1B̄1b1b̄1, lattice, x...; kwargs...)


# function _abba2aabb(x::GrassmannMPS, lattice::AbstractGrassmannLattice; trunc::TruncationScheme=DefaultKTruncation)
#     @assert length(x) == length(lattice)
#     x2 = deepcopy(x)
#     canonicalize!(x2, alg=Orthogonalize(trunc=trunc, normalize=false))
#     _abba2aabb_band!(x2, lattice, 0, trunc=trunc)
#     for j in lattice.k:-1:1
#         _abba2aabb_band!(x2, lattice, j, trunc=trunc)
#     end
#     return x2
# end

# function _abba2aabb_real(x::GrassmannMPS, lattice::AbstractGrassmannLattice; trunc::TruncationScheme=DefaultKTruncation)
#     @assert length(x) == length(lattice)
#     x2 = deepcopy(x)
#     canonicalize!(x2, alg=Orthogonalize(trunc=trunc, normalize=false))
#     _abba2aabb_band!(x2, lattice, 0, trunc=trunc)
#     for forward in (true, false)
#         # _abba2aabb_band_real!(x2, lattice, 0, trunc=trunc, forward=forward)
#         for j in lattice.k:-1:1
#             _abba2aabb_band_real!(x2, lattice, j, trunc=trunc, forward=forward)
#         end 
#     end
#     return x2
# end

# function _abba2aabb_band!(x::GrassmannMPS, lattice::AbstractGrassmannLattice, j::Int; trunc)
#     posa, posb = band_boundary(lattice, j)
#     _abba2aabb_band_util!(x, lattice, posa, posb, trunc=trunc)
# end

# function _abba2aabb_band_real!(x::GrassmannMPS, lattice::AbstractGrassmannLattice, j::Int; trunc, forward)
#     posa, posb = band_boundary(lattice, j, forward=forward)
#     _abba2aabb_band_util!(x, lattice, posa, posb, trunc=trunc)
# end

# function _abba2aabb_band_util!(x::GrassmannMPS, lattice::AbstractGrassmannLattice, posa::Int, posb::Int; trunc)
#     while posb >= posa + 1
#         for i in posb-1:-1:posa+1
#             easy_swap!(x, i, trunc=trunc)
#         end         
#         posa += 2
#     end 
# end

