module GTEMPO

# GrassmannTensor
export GrassmannTensorMap

# Grassmann MPS
export AbstractGTerm, GTerm, ExpGTerm, grassmannpspace
export AbstractGMPS, AbstractFiniteGMPS, GrassmannMPS, scaling, setscaling!, randomgmps, increase_bond!
export mult!, mult, DMRGMult1, DMRGMult2, DMRGMultAlgorithm
export GrassmannTransferMatrix

# Ordering of grasmann numbers
export GrassmannOrdering, ImagGrassmannOrdering, RealGrassmannOrdering, MixedGrassmannOrdering
export AbstractGrassmannLattice, ImagGrassmannLattice, RealGrassmannLattice, MixedGrassmannLattice, ContourIndex
export branches, matchindices, indexmappings, swapbandperm, swapband!, swapband
export OrderingStyle, ConjugationStyle, AdjacentConjugation, GeneralConjugation
export LayoutStyle, TimeLocalLayout, BandLocalLayout, BranchLocalLayout
export A1A1B1B1, AABB, A1B1B1A1, ABBA, A2A2A1A1B2B2B1B1
export A1A1B1B1a1a1b1b1, AABBaabb, A1A1a1a1B1B1b1b1, AAaaBBbb
export A1a1B1b1b1B1a1A1, AaBbbBaA, A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1, A1B1ā1b̄1A1B1a1b1
export A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2, A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2 #band local ordering
export A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2, AABB_AAaaBBbb, A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2 #mixedtime lattice
export A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2, AABB_aaAAbbBB
export changeordering, toadjacentordering

export vacuumstate, makestep, timesteps
export ImagGrassmannLattice1Order, RealGrassmannLattice1Order, RealGrassmannLattice2Order, GrassmannLattice, index, integrate
export IntegrationAlgorithm, ExactIntegrate, BMPSIntegrate, Zvalue

# correlation functions
export Cτ, Ct, Cm, branch, CorrelationMatrix, correlationfunction

# influence functional
export InfluenceFunctionalAlgorithm, PartialIF, TranslationInvariantIF
export partialinfluencefunctional
export influenceoperator, influenceoperatorexponential, differentialinfluencefunctional
export hybriddynamics, hybriddynamics!, hybriddynamicsstepper, hybriddynamicsstepper!

# GF and other observables
export gf, Gτ, parallel_Gτ, Gt, parallel_Gt, Gm, greater, lesser, contour_ordered_gf
export occupation, electriccurrent, electriccurrent_fast
export cached_gf, cached_Gτ, cached_Gt, cached_Gm, cached_greater, cached_lesser, cached_contour_ordered_gf
export cached_occupation, cached_electriccurrent, cached_electriccurrent_fast

# utilities for TEMPO
# exact models
export AbstractImpurityModel, SISB, SIDB, IRLM, SKIM, boundarycondition!, boundarycondition, boundarycondition_branching
export systhermalstate, systhermalstate!, sysdynamics, sysdynamics!, sysdynamicsstepper!, accsysdynamics, accsysdynamics_fast
export sysdynamics_forward!, sysdynamics_backward!, sysdynamics_imaginary!
export zoomin, zoomout
export ImpurityHamiltonian, tunneling, interaction, TunnelingTerm, InteractionTerm

using Base: @boundscheck, @propagate_inbounds
using Logging: @warn
using QuadGK, Permutations, Reexport, TupleTools, Strided, Statistics
using SphericalTensors: SphericalTensors, QR, SVD, LQ, AdjointTensorMap
const TK = SphericalTensors
using TensorOperations: TensorOperations, IndexTuple, Index2Tuple, linearize, Backend
const TO = TensorOperations
@reexport using DMRG, ImpurityModelBase
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

# correlation functions
include("correlations/correlations.jl")

# Feynman-Vernon influence functional as a multiplications of partial MPOs
include("influencefunctional/influencefunctional.jl")

# calculating observables and green's functions
include("gf/gf.jl")

# utility functions and models
include("models/models.jl")

end