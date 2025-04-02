include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


hybriddynamics(gmps::FockMPS, lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics(lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(vacuumstate(lattice), lattice, corr; kwargs...)


hybriddynamics_naive(gmps::FockMPS, lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics_naive(lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr; kwargs...)
