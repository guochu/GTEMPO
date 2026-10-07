include("svdmult.jl")
include("util.jl")
include("iterativemult.jl")

# wrapper

"""
    mult!(x::GrassmannMPS, y::GrassmannMPS, alg) -> (x, maxerr_or_info)

Multiplication of two GMPS x and y; the result is stored in x.
alg can be SVDCompression, DMRG1 or DMRG2

SVDCompression: standard SVD canonicalization, returns `(x, maxerr)`
DMRG1 and DMRG2: one- and two-site DMRG algorithm, returns `(x, info)` with `info`
the `ALSConvergenceInfo` provided by FiniteMPSAlgorithms

If one of x and y has very small bond dimension, then
SVDCompression is the method of choice,
when both x and y have large bond dimensions, then one
should use DMRGAlgorithm, and the perferred choice
is DMRG1
"""
mult(x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression) = mult(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
mult(x::GrassmannMPS, y::GrassmannMPS, alg::DMRGAlgorithm) = iterativemult(x, y, alg)

mult!(x::GrassmannMPS, y::GrassmannMPS, alg::SVDCompression) = mult!(x, y, trunc=alg.trunc, verbosity=alg.verbosity)
function mult!(x::GrassmannMPS, y::GrassmannMPS, alg::DMRGAlgorithm)
	z, info = iterativemult(x, y, alg)
	copy!(x.data, z.data)
	copy!(x.s, z.s)
	setscaling!(x, scaling(z))
	return x, info
end