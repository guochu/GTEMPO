include("svdmult.jl")
include("iterativemult.jl")

# wrapper
mult(x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression) = mult(x, y, trunc=alg.trunc)
mult(x::GrassmannMPS, y::GrassmannMPS, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)


const DefaultMultAlg = DMRGMult1(DefaultITruncation)