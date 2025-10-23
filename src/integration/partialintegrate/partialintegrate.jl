
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
partialintegrate(xs::GrassmannMPS...; cidx::Vector{Int}, trunc::TruncationScheme, verbosity::Int=0) = parint_mult(xs...; cidx=cidx, trunc=trunc, verbosity=verbosity)

function partialintegrate(lattice::AbstractGrassmannLattice, alg, xs::GrassmannMPS...; branchs::Tuple, bands::Tuple)
    unique(length.(xs)) == [length(lattice),] || throw(DimensionMismatch("unique($(length.(xs))) != [$(length(lattice)),]"))
    cidx = sort(unique([index(lattice, i; conj=false, band=band, branch=branch) for i in 0:lattice.k for band in bands for branch in branchs]))
    partialintegrate(alg, xs...; cidx=cidx)
end

