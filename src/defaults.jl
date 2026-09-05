#default settings
module Defaults
	const maxiter = 100 # for DMRG iteration
	const D = 100 # default bond dimension 
	const tolgauge = 1e-14 # for MPS truncation
	const tol = 1e-12 # for DMRG iteration
	const tollanczos = 1.0e-10 # for lanczos eigensolver
	const tolexp = 1.0e-8 # for local eigen in DMRG
	const verbosity = 1
end

const DefaultTruncation = truncdimcutoff(D=Defaults.D, ϵ=Defaults.tolgauge, add_back=0)

const DefaultIntegrationTruncation = truncdimcutoff(D=10000, ϵ=1.0e-12, add_back=0)
const DefaultITruncation = truncdimcutoff(D=200, ϵ=1.0e-10, add_back=0)
const DefaultKTruncation = truncdimcutoff(D=1000, ϵ=1.0e-10, add_back=0)
const DefaultMPOTruncation = truncdimcutoff(D=10000, ϵ=1.0e-10, add_back=0)

