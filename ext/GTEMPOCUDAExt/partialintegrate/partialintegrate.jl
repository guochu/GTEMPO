
include("svdmult.jl")
include("itermult.jl")


GTEMPO._ac_partialintegrate(alg::CuDMRGMultAlgorithm, xs::GrassmannMPS...; cidx::Vector{Int}) = cu_parint_iterativemult(xs...; cidx=cidx, alg=alg)
GTEMPO._ac_partialintegrate(alg::CuSVDCompression, xs::GrassmannMPS...; cidx::Vector{Int}) = cu_parint_mult(xs...; cidx=cidx, trunc=alg.trunc, verbosity=alg.verbosity)
