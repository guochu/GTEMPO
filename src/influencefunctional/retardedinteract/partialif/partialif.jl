include("util.jl")

# static mode, compute the whole IF and then measure
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")

retardedinteractdynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = retardedinteractdynamics!(vacuumstate(lattice), lattice, corr; kwargs...)
retardedinteractdynamics_naive(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = retardedinteractdynamics_naive!(vacuumstate(lattice), lattice, corr; kwargs...)
