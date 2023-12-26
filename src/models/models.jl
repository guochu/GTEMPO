abstract type AbstractImpurityModel end
sys_size(x::AbstractImpurityModel) = error("sys_size not implemented for model type $(typeof(x))")

hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityModel) = error("hybriddynamics not implemented for model $(typeof(model))")
hybriddynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; kwargs...) = hybriddynamics(vacuumstate(lattice), lattice, model; kwargs...)

sysdynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; kwargs...) = sysdynamics!(copy(gmps), lattice, model; kwargs...)
sysdynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; kwargs...) = error("sysdynamics! not implemented for model $(typeof(model))")
sysdynamicsstepper!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; kwargs...) = error("sysdynamicsstepper! not implemented for model $(typeof(model))")
sysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; kwargs...) = sysdynamics!(vacuumstate(lattice), lattice, model; kwargs...)



function sysdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityModel; 
						forward::Union{Nothing, Bool}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(forward)
		sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	else
		return forward ? sysdynamics_forward!(gmps, lattice, model; trunc=trunc) : sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	end
end 


sysdynamics_forward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityModel; kwargs...) = error("sysdynamics_forward! not implemented for model $(typeof(model))")
sysdynamics_backward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityModel; kwargs...) = error("sysdynamics_backward! not implemented for model $(typeof(model))")


systhermalstate(lattice::RealGrassmannLattice, model::AbstractImpurityModel; kwargs...) = systhermalstate!(vacuumstate(lattice), lattice, model; kwargs...)
"""
	systhermalstate!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityModel; kwargs...)
	Set the initial state of the impurity to be a locally thermal state
"""
function systhermalstate!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityModel; kwargs...)
	return systhermalstate_iterative!(gmps, lattice, model; kwargs...)
end 

include("boundary.jl")
include("correlationfunctions.jl")
include("hybriddynamics.jl")
include("sysdynamics.jl")
include("oneimpurity.jl")
include("irlm.jl")
include("skmodel.jl")
include("sysinitstate.jl")