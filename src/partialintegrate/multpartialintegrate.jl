

function check_contract_idx(L::Int, cidx::Vector{Int})
    @assert issorted(cidx) "contracted index should sorted"
    @assert all(isodd, cidx) "contracted index should be odd, which will contracted with adjacent one after it"
    @assert all(diff(cidx) .> 0) "cidx should be sorted, and no duplication"
    @assert (cidx[1] >= 1) && (cidx[end] < L)
end



_ac_partialintegrate(alg::DMRGMultAlgorithm, xs::Vararg{GrassmannMPS}; cidx::Vector{Int}) = parint_iterativemult(xs...; cidx=cidx, alg=alg)
_ac_partialintegrate(alg::SVDCompression, xs::Vararg{GrassmannMPS}; cidx::Vector{Int}) = parint_mult(xs...; cidx=cidx, trunc=alg.trunc, verbosity=alg.verbosity)
# _partialintegrate(xs::GrassmannMPS...; cidx::Vector{Int}, trunc::TruncationScheme, verbosity::Int=0) = parint_mult(xs...; cidx=cidx, trunc=trunc, verbosity=verbosity)

"""
	partialintegrate(lattice::AbstractGrassmannLattice, xs::GrassmannMPS...; branchs::Tuple, bands::Tuple, alg::DMRGAlgorithm=DefaultMultAlg)

!!! warning "Experimental interface"
    This is an experimental interface and is very likely to change
    significantly in future releases.

Integrate out all Grassmann variables on the given `branchs` and `bands`
(e.g. the bath branches `(:+, :)` on band 1), leaving only the system GVs.
`alg` (keyword, default `DefaultMultAlg`) selects the compression algorithm
(`SVDCompression`, `DMRG1`, ...).
"""
function partialintegrate(lattice::AbstractGrassmannLattice, xs::Vararg{GrassmannMPS}; branchs::Tuple, bands::Tuple, alg::DMRGAlgorithm=DefaultMultAlg)
    unique(length.(xs)) == [length(lattice),] || throw(DimensionMismatch("unique($(length.(xs))) != [$(length(lattice)),]"))
    cidx = sort(unique([index(lattice, i; conj=false, band=band, branch=branch) for i in 0:lattice.k for band in bands for branch in branchs]))
    _ac_partialintegrate(alg, xs...; cidx=cidx)
end

