push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/DMRG/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/InfiniteDMRG/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/GeneralHamiltonians/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/TEBD/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/FermionicTCMPS/src")

using Test, Random
using SphericalTensors, DMRG, FermionicTCMPS
const TK = SphericalTensors

# include("../src/includes.jl")

push!(LOAD_PATH, "../src")
using GTEMPO

Random.seed!(12354)

include("util.jl")


### tempo
# how to verify we get the correct correlation in an easy way?

const imag_grassmann_orderings = [A1A1B1B1(), A1B1B1A1(), A2A2A1A1B2B2B1B1()]
const real_grassmann_orderings = [A1A1a1a1B1B1b1b1(), A1a1B1b1b1B1a1A1(), A1B1ā1b̄1A1B1a1b1(), A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1(), A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2(), A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2()]
const imag_ac_grassmann_orderings = [A1A1B1B1(), A2A2A1A1B2B2B1B1()]
const real_ac_grassmann_orderings = [A1A1a1a1B1B1b1b1(), A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1(), A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2()]
const mixed_grassmann_orderings = [A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2(), A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2(), A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2()]
const mixed_ac_grassmann_orderings = [A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2(), A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2()]

include("tempo/grassmannlattice.jl")
include("tempo/grassmannmps.jl")
include("tempo/bmps_integrate.jl")

include("tempo/influencefunctional.jl")
include("tempo/influenceoperator.jl")
include("tempo/hybriddynamics.jl")

include("tempo/gf.jl")
include("tempo/cached_gf.jl")
include("tempo/bmps_gf.jl")

include("tempo/models.jl")
include("tempo/generalimpurity.jl")
include("tempo/observables.jl")

include("tempo/buildK.jl")


### tempo for interacting systems
include("Interacting/neq_tempo.jl")