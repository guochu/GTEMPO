include("svdmult.jl")
include("util.jl")
include("iterativemult.jl")

# wrapper

"""
    mult!(x::GrassmannMPS, y::GrassmannMPS, alg)

Multiplication of two GMPS x and y, and return the result 
alg can be SVDCompression, DMRG1 or DMRG2

SVDCompression: standard SVD canonicalization
DMRG1 and DMRG2: one- and two-site DMRG algorithm

If one of x and y has very small bond dimension, then 
SVDCompression is the method of choice,
when both x and y have large bond dimensions, then one 
should use DMRGMultAlgorithm, and the perferred choice 
is DMRG1
"""
mult(x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression) = mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
mult(x::GrassmannMPS, y::GrassmannMPS, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)

mult!(x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression) = mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
mult!(x::GrassmannMPS, y::GrassmannMPS, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)


const DefaultMultAlg = DMRG1(DefaultITruncation)