const DefaultIntegrationTruncation = truncdimcutoff(D=10000, ϵ=1.0e-12, add_back=0)
const DefaultITruncation = truncdimcutoff(D=200, ϵ=1.0e-7, add_back=0)
const DefaultKTruncation = truncdimcutoff(D=1000, ϵ=1.0e-10, add_back=0)
const DefaultMPOTruncation = truncdimcutoff(D=10000, ϵ=1.0e-10, add_back=0)
# a convenience function
DMRG.DMRG1(trunc::TruncationDimCutoff; maxiter::Int=5, tol::Float64=trunc.ϵ, verbosity::Int=0) = DMRG1(
			D=trunc.D, tolgauge=trunc.ϵ, maxiter=maxiter, tol=tol, verbosity=verbosity)
const DefaultMultAlg = DMRG1(DefaultITruncation)