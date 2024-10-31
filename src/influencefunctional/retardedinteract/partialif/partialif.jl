include("util.jl")

# static mode, compute the whole IF and then measure
include("imaginarytime.jl")
include("realtime.jl")


retardedinteractdynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = retardedinteractdynamics!(vacuumstate(lattice), lattice, corr; kwargs...)
