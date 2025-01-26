include("svdmult.jl")
include("iterativemult.jl")

# wrapper

"""
    mult!(x::GrassmannMPS, y::GrassmannMPS, alg)

Multiplication of two GMPS x and y, and return the result 
alg can be SVDCompression, DMRGMult1 or DMRGMult2

SVDCompression: standard SVD canonicalization
DMRGMult1 and DMRGMult2: one- and two-site DMRG algorithm

If one of x and y has very small bond dimension, then 
SVDCompression is the method of choice,
when both x and y have large bond dimensions, then one 
should use DMRGMultAlgorithm, and the perferred choice 
is DMRGMult1
"""
mult(x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression) = mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
mult(x::GrassmannMPS, y::GrassmannMPS, alg::DMRGMultAlgorithm) = iterativemult(x, y, alg)


const DefaultMultAlg = DMRGMult1(DefaultITruncation)