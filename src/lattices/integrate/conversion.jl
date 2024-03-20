# grassmann ordering conversion
# changeordering(::Type{A}, lattice::AbstractGrassmannLattice{A}, x::GrassmannMPS, ys::GrassmannMPS...; kwargs...) where A = x, ys...
# function changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, zs::GrassmannMPS...; kwargs...) where A
#     x1 = convert_ordering(A, lattice, x; kwargs...)
#     y1 = convert_ordering(A, lattice, y; kwargs...)
#     zs = map(z -> convert_ordering(A, lattice, z; kwargs...), zs)
#     return x1, y1, zs...
# end 
# changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS; kwargs...) where A = convert_ordering(A, lattice, x; kwargs...)

function changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, zs::GrassmannMPS...; 
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

# function convert_ordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation) where {A <: GrassmannOrdering} 
#     lattice2 = similar(lattice, ordering=A())
#     return permute(x, matchindices2(lattice2, lattice), trunc=trunc)
# end


# function convert_ordering(::Type{A1A1B1B1}, lattice::ImagGrassmannLattice{A1B1B1A1}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
#     return _abba2aabb(x, lattice, trunc=trunc)
# end

# function convert_ordering(::Type{A1A1a1a1B1B1b1b1}, lattice::RealGrassmannLattice{A1a1B1b1b1B1a1A1}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
#     return _abba2aabb(x, lattice, trunc=trunc)
# end

# function convert_ordering(::Type{A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, lattice::RealGrassmannLattice{A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
#     return _abba2aabb_real(x, lattice, trunc=trunc)
# end

# function convert_ordering(::Type{A1A1a1a1B1B1b1b1}, lattice::RealGrassmannLattice{A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
#     y = copy(x)
#     for j in 1:lattice.k
#         posa = index(lattice, j, band=1, conj=true, forward=true) + 1
#         posb = index(lattice, j, band=1, conj=false, forward=false) 
#         posb_end = index(lattice, j, band=lattice.bands, conj=false, forward=false) 
#         while (posb >= posa) && (posb <= posb_end)
#             for i in posb-1:-1:posa
#                 easy_swap!(y, i, trunc=trunc)
#             end
#             for i in posb:-1:posa+1
#                 easy_swap!(y, i, trunc=trunc)
#             end  
#             posa += 4
#             posb += 2         
#         end
#     end
#     return y
# end
# function convert_ordering(::Type{A1A1a1a1B1B1b1b1}, lattice::RealGrassmannLattice{A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, x::GrassmannMPS; kwargs...)
#     y = convert_ordering(A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2, lattice, x; kwargs...)
#     return convert_ordering(A1A1a1a1B1B1b1b1, similar(lattice, ordering=A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2()), y; kwargs...)
# end

"""
    toadjacentordering(lattice::AbstractGrassmannLattice, x::GrassmannMPS...; kwargs...)

Convert the lattice into a "closest" adjacent ordering
This function a simple wrapper of changeordering, by specifying a particular ordering 
"""
toadjacentordering(lattice::AbstractGrassmannLattice, x::GrassmannMPS...; kwargs...) = error("toadjacentordering not implemented for lattice type $(typeof(lattice))")
toadjacentordering(lattice::ImagGrassmannLattice{<:A1B1B1A1}, x::GrassmannMPS...; kwargs...) = changeordering(A1A1B1B1, lattice, x...; kwargs...)
# toadjacentordering(lattice::RealGrassmannLattice{<:A1a1B1b1b1B1a1A1}, x::GrassmannMPS...; kwargs...) = changeordering(A1A1a1a1B1B1b1b1, lattice, x...; kwargs...)
# toadjacentordering(lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, x::GrassmannMPS...; kwargs...) = changeordering(
#                     A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2, lattice, x...; kwargs...)
function toadjacentordering(lattice::RealGrassmannLattice, x::GrassmannMPS...; kwargs...)
    if LayoutStyle(lattice) isa TimeLocalLayout
        return changeordering(A1A1a1a1B1B1b1b1, lattice, x...; kwargs...)
    else
        return changeordering(A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2, lattice, x...; kwargs...)
    end
end 

toadjacentordering(lattice::MixedGrassmannLattice{<:A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, x::GrassmannMPS...; kwargs...) = changeordering(
                    A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2, lattice, x...; kwargs...)

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

