using Test, Random
using GTEMPO
using Z2Tensors
using ImpurityModelBase

Random.seed!(12354)

# shared utilities: spectra, Grassmann orderings, ED reference builders

using GTEMPO: A1B1B̄1Ā1, Ā2A1B̄2B1,
          A1B1ā1b̄1Ā1B̄1a1b1, A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2,
          Ā2A1ā1a2B̄2B1b̄1b̄2, A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2,
          Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2, A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1,
          FockMatrix, fock_propagator,
          fock_thermalstate, boundarycondition_branching

include("util.jl")

# Part 1: basic API tests (fast sanity checks of the public interface)
include("api/runtests.jl")

# Part 2: functional tests against exact references (ED / analytic)
include("normalbath/runtests.jl")

# BCS / electron-phonon code paths (few-mode & continuous bath)
include("bcs/runtests.jl")
include("electronphonon/runtests.jl")

# Part 3: multi-orbital impurity models vs exact diagonalization
include("multiorbital/runtests.jl")
