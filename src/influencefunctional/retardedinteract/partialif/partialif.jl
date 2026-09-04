include("util.jl")

# static mode, compute the whole IF and then measure
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")

"""
	retardedinteractdynamics_naive(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...)

Build the retarded interaction e^{ΣᵢⱼΔᵢⱼnᵢnⱼ} as a GMPS using the PartialIF algorithm
(the TTI-IF algorithm for this case is to be developed)
corr: bosonic hybridization function calculated using QuAPI
"""
retardedinteractdynamics_naive(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = retardedinteractdynamics_naive!(vacuumstate(lattice), lattice, corr; kwargs...)
