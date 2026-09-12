module GTEMPO

# auxiliary 
export truncdimcutoff

# GrassmannTensor
export scalartype
export @grassmann
using Z2Tensors: NoTruncation
export NoTruncation

# Grassmann MPS
export AbstractGTerm, GTerm, ExpGTerm, z2space
export AbstractGMPS, AbstractFiniteGMPS, GrassmannMPS, SparseGMPS, togmps, scaling, setscaling!, randomgmps, increase_bond!
export iscanonical, isleftcanonical, isrightcanonical
export mult!, mult, DMRG1, DMRG2
export GrassmannTransferMatrix
export randomfockmps

# definition of lattice and Ordering of grassmann numbers
export GrassmannOrdering, ImagGrassmannOrdering, RealGrassmannOrdering, MixedGrassmannOrdering
export AbstractGrassmannLattice, ImagGrassmannLattice, RealGrassmannLattice, MixedGrassmannLattice, ContourIndex
export branches, matchindices, indexmappings, swapbandperm, swapband!, swapband, fillband
export OrderingStyle, ConjugationStyle, AdjacentConjugation, GeneralConjugation
export LayoutStyle, TimeLocalLayout, BranchLocalLayout, GeneralLayout
# export TimeOrderingStyle, ImaginaryTimeOrderingStyle, RealTimeOrderingStyle, TimeAscending, TimeDscending
export A1Ā1B1B̄1
export A1Ā1B1B̄1a1ā1b1b̄1
export A1Ā1a1ā1B1B̄1b1b̄1
export A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2
export A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2
export A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2
export vacuumstate, makestep, timesteps
export ImagGrassmannLattice1Order, RealGrassmannLattice1Order, RealGrassmannLattice2Order, GrassmannLattice, index

# integration of GMPSs
export integrate, integrateband, integratebands, partialintegrate, multintegrateband
export IntegrationAlgorithm, ExactIntegrate, BMPSIntegrate, Zvalue
export changeordering, toadjacentordering
export environments

# correlation functions
export branch, correlationfunction

# influence functional
export InfluenceFunctionalAlgorithm, PartialIF, XTRGIF, ExactTTIIF, TDVPIF, partialif_hybrid, partialif_hybrid_naive
export influenceoperators, influenceoperatorsteppers, influenceoperatorstepper
export hybriddynamics, hybriddynamics!, hybriddynamics_naive, hybriddynamics_naive!, hybriddynamicsstepper, hybriddynamicsstepper!
export retardedinteractdynamics_naive, retardedinteractdynamics_naive!

# GF and other observables
export gf, greater, lesser, contour_ordered_gf
export occupation, electriccurrent, electriccurrent_fast, heatcorrelationfunction, heatcurrent_fast
export cached_gf, cached_greater, cached_lesser, cached_contour_ordered_gf
export cached_occupation, cached_electriccurrent, cached_electriccurrent_fast, cached_heatcurrent_fast
export cached_gf_fast
export cached_greater_fast, cached_lesser_fast
export nn, cached_nn, insert_n!, insert_n, nn2, cached_nn2

# connections of Grassmann variables
export bulkconnection!, bulkconnection
export boundarycondition!, boundarycondition

# utilities for TEMPO
# impurity model Hamilltonians
export AbstractImpurityHamiltonian, ConstImpurityHamiltonian, AbstractTdImpurityHamiltonian, num_bands
export AndersonIM, ToulouseIM, IRLM, KanamoriIM
export sysdynamics, sysdynamics_fast, sysdynamicsstepper!
export sysdynamics_imaginary!, sysdynamics_forward!, sysdynamics_backward!
export sysinitialstate, sysinitialstate!, systhermalstate, systhermalstate!
export ImpurityHamiltonian, QuenchedImpurityHamiltonian, TdImpurityHamiltonian, TdImpurityOp
export baresysdynamics, baresysdynamics_fast
export baresysdynamics_imaginary!, baresysdynamics_forward!, baresysdynamics_backward!
export sysdynamics_deprecated, sysdynamics_deprecated!
export sysdynamics_imaginary_deprecated!, sysdynamics_forward_deprecated!, sysdynamics_backward_deprecated!
export baresysdynamics_deprecated, baresysdynamics_deprecated!


# electron-phonon interaction
export AbstractNTerm, ExpNTerm
export FockMPS
export FockOrdering, ImagFockOrdering, RealFockOrdering, MixedFockOrdering, similargrassmannordering
export M1N1, MN, M1m1N1n1, MmNn, M1N1_M1m1N1n1M2m2N2n2, MN_MmNn, M1N1_m1M1n1N1m2M2n2N2
export AbstractFockLattice, FockLattice, ImagFockLattice, similargrassmannlattice
export RealFockLattice, MixedFockLattice
export reweighting!, reweighting


export MPO, PartialMPO, AbstractMPO, MPOHamiltonian, SchurMPOTensor, apply!
export MPSTensor, MPSBondTensor, MPOTensor, SiteOperator
export canonicalize, canonicalize!, physical_spaces, environments, expectationvalue, positions, physical_space
export ophysical_space, iphysical_space
export svectors_uninitialized, unset_svectors!, mpotensortype, mpstensortype
export SVDCompression
export DMRGAlgorithm
export timeevompo, WI, WII, ComplexStepper, FirstOrderStepper, complex_stepper
export bond_dimension, bond_dimensions, distance, space_l, space_r, l_LL, r_RR, bondtensortype
export Orthogonalize
export AbstractPronyExpansion, OverDeterminedProny, DeterminedProny, MatrixPencil, LeastSquareProny



using Base: @boundscheck, @propagate_inbounds
using Logging: @warn
using Reexport, TupleTools, Strided
using Z2Tensors
using Z2Tensors: Z2Tensors, QR, SVD, LQ, AdjointTensorMap, NoTruncation, TruncationDimCutoff
const TK = Z2Tensors
using TensorOperations: TensorOperations, IndexTuple, Index2Tuple, linearize, AbstractBackend # for Grassmann Tensors
const TO = TensorOperations
# @reexport using DMRG, ImpurityModelBase, QuAPI
@reexport using ImpurityModelBase, QuAPI
import QuAPI: branch, index
# using DMRG: TimeEvoMPOAlgorithm






using ExpExp
using KrylovKit: Arnoldi, exponentiate
using LinearAlgebra: LinearAlgebra, Symmetric, eigen, qr, pinv, eigvals, Diagonal, diagm


include("auxiliary/linalg.jl")
include("auxiliary/mpstensors.jl")
include("algorithms.jl")
include("defaults.jl") # default constants (DefaultMultAlg requires DMRG1 from algorithms.jl)

include("mpo/mpo.jl")


# # TEMPO algorithm

# GrassmannTensor
include("grassmanntensor/grassmanntensor.jl")
include("grassmanntensor/tensoroperations.jl")
include("grassmanntensor/grassmannmacro.jl")

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
include("grassmannmps/sparsegmps.jl")

# Grassmann lattice
include("lattices/lattices.jl")

# integration
include("integration/integration.jl")

# partial integration
include("partialintegrate/partialintegrate.jl")

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
include("sysdynamics/sysdynamics.jl")

# electron phonon interactions
include("electronphonon/electronphonon.jl")
end