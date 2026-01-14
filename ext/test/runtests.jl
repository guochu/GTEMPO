ENV["CUDA_VISIBLE_DEVICES"] = 6

push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/DMRG/src")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/InfiniteDMRG/src")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/GeneralHamiltonians/src")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/TEBD/src")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/FermionicTCMPS/src")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/ImpurityModelBase/src")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/QuAPI/src")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/Z2Tensors/src")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/Z2Tensors/ext/Z2TensorsCUDAExt/")
push!(LOAD_PATH, dirname(dirname(dirname(Base.@__DIR__))) * "/GTEMPO/src")

using Test, Random
using Z2Tensors
const TK = Z2Tensors
using GTEMPO
                    

include("../GTEMPOCUDAExt/includes.jl")

Random.seed!(12354)

# include("util.jl")
const imag_grassmann_orderings = [A1Ā1B1B̄1(), A1B1B̄1Ā1(), A2Ā2A1Ā1B2B̄2B1B̄1()]
const real_grassmann_orderings = [A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A1B1ā1b̄1Ā1B̄1a1b1(), A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1(), A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2()]
const imag_ac_grassmann_orderings = [A1Ā1B1B̄1(), A2Ā2A1Ā1B2B̄2B1B̄1()]
const real_ac_grassmann_orderings = [A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(), A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2()]
const mixed_grassmann_orderings = [A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2(), A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2(), A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2()]
const mixed_ac_grassmann_orderings = [A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2(), A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2()]


include("grassmannmps.jl")
include("partialintegrate.jl")
include("cached_gf_fast2.jl")
