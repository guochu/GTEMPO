using Base: @boundscheck, @propagate_inbounds
using Logging: @warn
using Permutations, Reexport, TupleTools, Strided, Statistics, TensorKit
using TensorKit: TensorKit, QR, SVD, LQ, AdjointTensorMap, NoTruncation
const TK = TensorKit
using TensorOperations: TensorOperations, IndexTuple, Index2Tuple, linearize, AbstractBackend
const TO = TensorOperations
@reexport using DMRG, ImpurityModelBase, QuAPI
import QuAPI: branch, index
using DMRG: TimeEvoMPOAlgorithm


# # TEMPO algorithm

# GrassmannTensor
include("grassmanntensor/grassmanntensor.jl")
include("grassmanntensor/linalg.jl")
include("grassmanntensor/tensoroperations.jl")

# default constants
include("defaults.jl")

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

# Grassmann lattice and integration
include("lattices/lattices.jl")

# # correlation functions
include("correlationfunction.jl")

# Feynman-Vernon influence functional as a multiplications of partial MPOs
include("influencefunctional/influencefunctional.jl")

# calculating observables and green's functions
include("gf/gf.jl")

# utility functions and models
include("models/models.jl")