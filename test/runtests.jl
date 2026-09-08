using Test, Random
using GTEMPO
using Z2Tensors
using ImpurityModelBase

Random.seed!(12354)

# shared utilities: spectra, Grassmann orderings, ED reference builders
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
