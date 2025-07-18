module GTEMPO

# GrassmannTensor
export GrassmannTensorMap

# Grassmann MPS
export AbstractGTerm, GTerm, ExpGTerm, grassmannpspace
export AbstractGMPS, AbstractFiniteGMPS, GrassmannMPS, scaling, setscaling!, randomgmps, increase_bond!
export mult!, mult, DMRGMult1, DMRGMult2, DMRGMultAlgorithm
export GrassmannTransferMatrix

# definition of lattice and Ordering of grassmann numbers
export GrassmannOrdering, ImagGrassmannOrdering, RealGrassmannOrdering, MixedGrassmannOrdering
export AbstractGrassmannLattice, ImagGrassmannLattice, RealGrassmannLattice, MixedGrassmannLattice, ContourIndex
export branches, matchindices, indexmappings, swapbandperm, swapband!, swapband, fillband
export OrderingStyle, ConjugationStyle, AdjacentConjugation, GeneralConjugation
export LayoutStyle, TimeLocalLayout, BandLocalLayout, BranchLocalLayout
# export TimeOrderingStyle, ImaginaryTimeOrderingStyle, RealTimeOrderingStyle, TimeAscending, TimeDscending
export A1Ā1B1B̄1, AĀBB̄, A1B1B̄1Ā1, ABB̄Ā, A2Ā2A1Ā1B2B̄2B1B̄1
export A1Ā1B1B̄1a1ā1b1b̄1, AĀBB̄aābb̄, A1Ā1a1ā1B1B̄1b1b̄1, AĀaāBB̄bb̄
export A1Ā1B1B̄1b̄1B̄1ā1Ā1, AaBbb̄B̄āĀ, A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1, ABāb̄ĀB̄ab, A1B1ā1b̄1Ā1B̄1a1b1
export A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2, A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2 #band local ordering
export A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2, AĀBB̄_AĀaāBB̄bb̄, A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2 #mixedtime lattice
export A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2, AĀBB̄_aāAĀbb̄BB̄
export vacuumstate, makestep, timesteps
export ImagGrassmannLattice1Order, RealGrassmannLattice1Order, RealGrassmannLattice2Order, GrassmannLattice, index

# integration of GMPSs
export integrate, integrateband
export IntegrationAlgorithm, ExactIntegrate, BMPSIntegrate, Zvalue
export changeordering, toadjacentordering

# correlation functions
export branch, correlationfunction

# influence functional
export InfluenceFunctionalAlgorithm, PartialIF, TranslationInvariantIF, ExactTranslationInvariantIF, partialif_hybrid, partialif_hybrid_naive
export influenceoperator, influenceoperatorexponential, differentialinfluencefunctional
export hybriddynamics, hybriddynamics!, hybriddynamics_naive, hybriddynamics_naive!, hybriddynamicsstepper, hybriddynamicsstepper!
export retardedinteractdynamics, retardedinteractdynamics!, retardedinteractdynamics_naive, retardedinteractdynamics_naive!

# GF and other observables
export gf, Gτ, parallel_Gτ, Gt, parallel_Gt, Gm, greater, lesser, contour_ordered_gf
export occupation, occupation2, electriccurrent, electriccurrent_fast, heatcurrent_fast
export cached_gf, cached_Gτ, cached_Gt, cached_Gm, cached_greater, cached_lesser, cached_contour_ordered_gf
export cached_occupation, cached_electriccurrent, cached_electriccurrent_fast, cached_heatcurrent_fast
export cached_gf_fast, cached_Gτ_fast, cached_Gt_fast, cached_Gm_fast
export cached_greater_fast, cached_lesser_fast
export nn, cached_nn, insert_n!, insert_n, nn2, cached_nn2

# connections of Grassmann variables
export bulkconnection!, bulkconnection,
export boundarycondition!, boundarycondition, boundarycondition_branching

# utilities for TEMPO
# impurity model Hamilltonians
export AbstractImpurityHamiltonian, AndersonIM, IRLM, KanamoriIM
export systhermalstate, systhermalstate!, sysdynamics, sysdynamics!, sysdynamicsstepper!, accsysdynamics, accsysdynamics_fast
# export sysdynamics_forward!, sysdynamics_backward!, sysdynamics_imaginary!
export zoomin, zoomout
export ImpurityHamiltonian, tunneling, interaction, TunnelingTerm, InteractionTerm, AbstractFTerm
export baresysdynamics!, baresysdynamics


# electron-phonon interaction
export AbstractNTerm, ExpNTerm
export FockMPS
export FockOrdering, ImagFockOrdering, RealFockOrdering, MixedFockOrdering, similargrassmannordering
export M1N1, MN, M1m1N1n1, MmNn, M1N1_M1m1N1n1M2m2N2n2, MN_MmNn
export AbstractFockLattice, FockLattice, ImagFockLattice, similargrassmannlattice
export RealFockLattice, MixedFockLattice
export reweighting!, reweighting

using Base: @boundscheck, @propagate_inbounds
using Logging: @warn
using Permutations, Reexport, TupleTools, Strided, Statistics, TensorKit
using TensorKit: TensorKit, QR, SVD, LQ, AdjointTensorMap, NoTruncation
const TK = TensorKit
using TensorOperations: TensorOperations, IndexTuple, Index2Tuple, linearize, AbstractBackend # for Grassmann Tensors
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

# Grassmann lattice
include("lattices/lattices.jl")

# integration
include("integration/integration.jl")

# # correlation functions
include("correlationfunction.jl")

# Feynman-Vernon influence functional as a multiplications of partial MPOs
include("influencefunctional/influencefunctional.jl")

# calculating observables and green's functions
include("observables/observables.jl")

# grassmann variables connections
include("gvconnections/gvconnections.jl")

# utility functions and models
include("models/models.jl")

# electron phonon interactions
include("electronphonon/electronphonon.jl")
end