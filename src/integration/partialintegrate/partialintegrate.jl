
include("utils.jl")
include("svdmult.jl")
include("itermult.jl")


function chack_contract_idx(L::Int, cidx::Vector{Int})
    @assert issorted(cidx) "contracted index should sorted"
    @assert all(isodd, cidx) "contracted index should be odd, which will contracted with adjacent one after it"
    @assert all(diff(cidx) .> 0) "cidx should be sorted, and no duplication"
    @assert (cidx[1] >= 1) && (cidx[end] < L)
end



partialintegrate(alg::DMRGMultAlgorithm, xs::GrassmannMPS...; cidx::Vector{Int}) = parint_iterativemult(xs...; cidx=cidx, alg=alg)
partialintegrate(alg::SVDCompression, xs::GrassmannMPS...; cidx::Vector{Int}) = parint_mult(xs...; cidx=cidx, trunc=alg.trunc, verbosity=alg.verbosity)
function partialintegrate(lattice::AbstractGrassmannLattice, bands::Tuple, alg, xs::GrassmannMPS...)
    unique(length.(xs)) == [length(lattice),] || throw(DimensionMismatch())
    cidx = sort([index(lattice, i; conj=false, band=band) for i in 0:lattice.k for band in bands])
    partialintegrate(alg, xs...; cidx=cidx)
end

partialintegrate(xs::GrassmannMPS...; cidx::Vector{Int}, trunc::TruncationScheme, verbosity::Int=0) = parint_mult(xs...; cidx=cidx, trunc=trunc, verbosity=verbosity)
partialintegrate(lattice::AbstractGrassmannLattice, bands::Tuple, xs::GrassmannMPS...; trunc::TruncationScheme, verbosity::Int=0) = partialintegrate(lattice, bands, SVDCompression(trunc; verbosity=verbosity), xs...)

