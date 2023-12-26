# grassmann ordering conversion
convert_ordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS; kwargs...) where A = error("convert_ordering not implemented for source $(typeof(lattice.ordering)) to target $(A)")
changeordering(::Type{A}, lattice::AbstractGrassmannLattice{A}, x::GrassmannMPS, ys::GrassmannMPS...; kwargs...) where A = x, ys...
function changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS, y::GrassmannMPS, zs::GrassmannMPS...; kwargs...) where A
    x1 = convert_ordering(A, lattice, x; kwargs...)
    y1 = convert_ordering(A, lattice, y; kwargs...)
    zs = map(z -> convert_ordering(A, lattice, z; kwargs...), zs)
    return x1, y1, zs...
end 
changeordering(::Type{A}, lattice::AbstractGrassmannLattice, x::GrassmannMPS; kwargs...) where A = convert_ordering(A, lattice, x; kwargs...)

function convert_ordering(::Type{A1A1B1B1}, lattice::ImagGrassmannLattice{A1B1B1A1}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
    return _abba2aabb(x, lattice, trunc=trunc)
end

function convert_ordering(::Type{A1A1a1a1B1B1b1b1}, lattice::RealGrassmannLattice{A1a1B1b1b1B1a1A1}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
    return _abba2aabb(x, lattice, trunc=trunc)
end

function convert_ordering(::Type{A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, lattice::RealGrassmannLattice{A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
    return _abba2aabb_real(x, lattice, trunc=trunc)
end

function convert_ordering(::Type{A1A1a1a1B1B1b1b1}, lattice::RealGrassmannLattice{A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, x::GrassmannMPS; trunc::TruncationScheme=DefaultKTruncation)
    y = copy(x)
    for j in 1:lattice.k
        posa = index(lattice, j, band=1, conj=true, forward=true) + 1
        posb = index(lattice, j, band=1, conj=false, forward=false) 
        posb_end = index(lattice, j, band=lattice.bands, conj=false, forward=false) 
        while (posb >= posa) && (posb <= posb_end)
            for i in posb-1:-1:posa
                easy_swap!(y, i, trunc=trunc)
            end
            for i in posb:-1:posa+1
                easy_swap!(y, i, trunc=trunc)
            end  
            posa += 4
            posb += 2         
        end
    end
    return y
end
function convert_ordering(::Type{A1A1a1a1B1B1b1b1}, lattice::RealGrassmannLattice{A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, x::GrassmannMPS; kwargs...)
    y = convert_ordering(A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2, lattice, x; kwargs...)
    return convert_ordering(A1A1a1a1B1B1b1b1, similar(lattice, ordering=A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2()), y; kwargs...)
end

toadjacentordering(lattice::AbstractGrassmannLattice, x::GrassmannMPS...; kwargs...) = error("toadjacentordering not implemented for lattice type $(typeof(lattice))")

function toadjacentordering(lattice::ImagGrassmannLattice{<:A1B1B1A1}, x::GrassmannMPS...; trunc::TruncationScheme=DefaultKTruncation)
    lattice2 = similar(lattice, ordering=A1A1B1B1())
    return lattice2, map(y -> _abba2aabb(y, lattice, trunc=trunc), x)...
end
function toadjacentordering(lattice::RealGrassmannLattice{<:A1a1B1b1b1B1a1A1}, x::GrassmannMPS...; trunc::TruncationScheme=DefaultKTruncation)
    lattice2 = similar(lattice, ordering=A1A1a1a1B1B1b1b1())
    return lattice2, map(y -> _abba2aabb(y, lattice, trunc=trunc), x)...
end
function toadjacentordering(lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, x::GrassmannMPS...; trunc::TruncationScheme=DefaultKTruncation)
    lattice2 = similar(lattice, ordering=A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2())
    return lattice2, map(y -> _abba2aabb_real(y, lattice, trunc=trunc), x)...
end

function _abba2aabb(x::GrassmannMPS, lattice::AbstractGrassmannLattice; trunc::TruncationScheme=DefaultKTruncation)
    @assert length(x) == length(lattice)
    x2 = copy(x)
    _abba2aabb_band!(x2, lattice, 0, trunc=trunc)
    for j in lattice.k:-1:1
        _abba2aabb_band!(x2, lattice, j, trunc=trunc)
    end
    return x2
end

function _abba2aabb_real(x::GrassmannMPS, lattice::AbstractGrassmannLattice; trunc::TruncationScheme=DefaultKTruncation)
    @assert length(x) == length(lattice)
    x2 = copy(x)
    _abba2aabb_band!(x2, lattice, 0, trunc=trunc)
    for forward in (true, false)
        # _abba2aabb_band_real!(x2, lattice, 0, trunc=trunc, forward=forward)
        for j in lattice.k:-1:1
            _abba2aabb_band_real!(x2, lattice, j, trunc=trunc, forward=forward)
        end 
    end
    return x2
end

function _abba2aabb_band!(x::GrassmannMPS, lattice::AbstractGrassmannLattice, j::Int; trunc)
    posa, posb = band_boundary(lattice, j)
    _abba2aabb_band_util!(x, lattice, posa, posb, trunc=trunc)
end

function _abba2aabb_band_real!(x::GrassmannMPS, lattice::AbstractGrassmannLattice, j::Int; trunc, forward)
    posa, posb = band_boundary(lattice, j, forward=forward)
    _abba2aabb_band_util!(x, lattice, posa, posb, trunc=trunc)
end

# function _abba2aabb_band!(x::GrassmannMPS, lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, j::Int; forward::Bool, trunc)
#     posa, posb = band_boundary(lattice, j, forward=forward)
#     abba2aabb_band_util!(x, lattice, posa, posb, trunc=trunc)
# end

function _abba2aabb_band_util!(x::GrassmannMPS, lattice::AbstractGrassmannLattice, posa::Int, posb::Int; trunc)
    while posb >= posa + 1
        for i in posb-1:-1:posa+1
            easy_swap!(x, i, trunc=trunc)
        end         
        posa += 2
    end 
end


function band_boundary(lattice::ImagGrassmannLattice{<:A1B1B1A1}, j::Int)
    posa = index(lattice, j, conj=false, band=1)
    posb = index(lattice, j, conj=true, band=1)     
    return posa, posb
end
function band_boundary(lattice::ImagGrassmannLattice{<:A1A1B1B1}, j::Int)
    posa = index(lattice, j, conj=false, band=1)
    posb = index(lattice, j, conj=true, band=lattice.bands)     
    return posa, posb
end

function band_boundary(lattice::RealGrassmannLattice{<:A1a1B1b1b1B1a1A1}, j::Int)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=1)     
    else
        posa = index(lattice, j, conj=false, forward=true, band=1)
        posb = index(lattice, j, conj=true, forward=true, band=1)
    end
    return posa, posb
end
function band_boundary(lattice::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, j::Int)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=lattice.bands)     
    else
        posa = index(lattice, j, conj=false, forward=true, band=1)
        posb = index(lattice, j, conj=true, forward=false, band=lattice.bands)
    end
    return posa, posb
end
function band_boundary(lattice::RealGrassmannLattice{<:A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2}, j::Int; forward::Union{Nothing, Bool}=nothing)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=1)   
    else
        @assert isa(forward, Bool)
        posa = index(lattice, j, conj=false, band=1, forward=forward)
        posb = index(lattice, j, conj=true, band=1, forward=forward) 
    end
    return posa, posb
end
function band_boundary(lattice::RealGrassmannLattice{<:A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2}, j::Int; forward::Union{Nothing, Bool}=nothing)
    if j == 0
        posa = index(lattice, j, conj=false, band=1)
        posb = index(lattice, j, conj=true, band=lattice.bands)   
    else
        @assert isa(forward, Bool)
        posa = index(lattice, j, conj=false, band=1, forward=forward)
        posb = index(lattice, j, conj=true, band=lattice.bands, forward=forward) 
    end
    return posa, posb
end
