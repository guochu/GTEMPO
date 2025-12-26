using Base: @boundscheck, @propagate_inbounds
using Logging: @warn
using Permutations, Reexport, TupleTools, Strided, Statistics
using Z2Tensors
using Z2Tensors: Z2Tensors, QR, SVD, LQ, AdjointTensorMap, NoTruncation, TruncationDimCutoff
const TK = Z2Tensors
using TensorOperations: TensorOperations, IndexTuple, Index2Tuple, linearize, AbstractBackend # for Grassmann Tensors
const TO = TensorOperations
# @reexport using DMRG, ImpurityModelBase, QuAPI
@reexport using ImpurityModelBase, QuAPI
import QuAPI: branch, index
# using DMRG: TimeEvoMPOAlgorithm










using Parameters, Polynomials, KrylovKit, LsqFit
using LinearAlgebra: LinearAlgebra, Symmetric, eigen, qr, pinv, eigvals, Diagonal

include("auxiliary/CachedVectors.jl")
include("auxiliary/defaults.jl") # default constants
include("auxiliary/linalg.jl")
include("auxiliary/mpstensors.jl")
include("auxiliary/orth.jl")
include("auxiliary/mpsalgs.jl")

include("mpo/mpo.jl")





# # TEMPO algorithm

# GrassmannTensor
include("grassmanntensor/grassmanntensor.jl")
include("grassmanntensor/linalg.jl")
include("grassmanntensor/tensoroperations.jl")

# Grassmann MPS operations
include("grassmannmps/util.jl")
include("grassmannmps/space.jl")
include("grassmannmps/grassmannterms.jl")
include("grassmannmps/abstractgmps.jl")
include("grassmannmps/grassmannmps.jl")
include("grassmannmps/orth.jl")
include("grassmannmps/linalg.jl")
include("grassmannmps/transfer.jl")
include("grassmannmps/mult/mult.jl")

# Grassmann lattice
include("lattices/lattices.jl")

# integration
include("integration/integration.jl")

# # correlation functions
include("correlationfunction.jl")

# Feynman-Vernon influence functional as a multiplications of partial MPOs
include("influencefunctional/influencefunctional.jl")

include("bcsinfluencefunctional/bcsinfluencefunctional.jl")

# calculating observables and green's functions
include("observables/observables.jl")

# grassmann variables connections
include("gvconnections/gvconnections.jl")

# utility functions and models
include("models/models.jl")

# electron phonon interactions
include("electronphonon/electronphonon.jl")