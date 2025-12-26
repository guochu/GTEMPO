push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/ImpurityModelBase/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/QuAPI/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/Z2TensorKit/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/Z2TensorKit/ext/Z2TensorKitCUDAExt/")

using Test, Random
using Z2TensorKit
const TK = Z2TensorKit
                    

include("../src/includes.jl")

# push!(LOAD_PATH, "../src")
# using GTEMPO

Random.seed!(12354)

include("util.jl")


# include("DMRG/tests.jl")




### tempo
# how to verify we get the correct correlation in an easy way?

const imag_grassmann_orderings = [A1Ā1B1B̄1(), A1B1B̄1Ā1(), A2Ā2A1Ā1B2B̄2B1B̄1()]
const real_grassmann_orderings = [A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A1B1ā1b̄1Ā1B̄1a1b1(), A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1(), A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2()]
const imag_ac_grassmann_orderings = [A1Ā1B1B̄1(), A2Ā2A1Ā1B2B̄2B1B̄1()]
const real_ac_grassmann_orderings = [A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(), A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2()]
const mixed_grassmann_orderings = [A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2(), A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2(), A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2()]
const mixed_ac_grassmann_orderings = [A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2(), A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2()]

include("grassmanntensor.jl")
include("tempo/grassmannlattice.jl")
include("tempo/grassmannmps.jl")
include("tempo/bmps_integrate.jl")

include("tempo/partialif_hybrid.jl")
include("tempo/retardedinteract.jl")
include("tempo/influenceoperator.jl")
include("tempo/hybriddynamics.jl")

include("tempo/gf.jl")
include("tempo/cached_gf.jl")
include("tempo/cached_nn.jl")
include("tempo/cached_gf_fast.jl")
include("tempo/cached_gf_fast2.jl")
include("tempo/bmps_gf.jl")
include("tempo/observables.jl")

include("tempo/swapband.jl")
include("tempo/fillband.jl")
include("tempo/integrateband.jl")
include("tempo/integrateband2.jl")
include("tempo/partialintegrate.jl")

include("tempo/models.jl")
include("tempo/independentbosons.jl")
include("tempo/irlm.jl")
include("tempo/generalimpurity.jl")
include("tempo/buildK.jl")

include("tempo/bcs/bcs.jl")

### tempo for interacting systems
include("interacting/neq_tempo.jl")

### electron-phonon interaction
include("electronphonon/focklattice.jl")
include("electronphonon/fockmps.jl")
include("electronphonon/retardedinteract.jl")
include("electronphonon/independentbosons.jl")


