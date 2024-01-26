include("svdmult.jl")
include("iterativemult.jl")

# wrapper
mult(x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression) = mult(x, y, trunc=alg.trunc)
mult(x::GrassmannMPS, y::GrassmannMPS, alg::Union{DMRG1, DMRG2}) = iterativemult(x, y, alg)