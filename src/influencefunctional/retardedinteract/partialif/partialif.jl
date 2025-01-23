include("util.jl")

# static mode, compute the whole IF and then measure
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")

"""
	hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg; band::Int=1) 

Build the retarded interaction e^{ΣᵢⱼΔᵢⱼnᵢnⱼ} as a GMPS using the PartialIF algorithm 
(the TTI-IF algorithm for this case is to be developed)
corr: bosonic hybridization function calculated using QuAPI
"""
retardedinteractdynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = retardedinteractdynamics!(vacuumstate(lattice), lattice, corr; kwargs...)
retardedinteractdynamics_naive(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = retardedinteractdynamics_naive!(vacuumstate(lattice), lattice, corr; kwargs...)
