# interface
# only the impurity Hamiltonian, no bath
abstract type AbstractImpurityHamiltonian end
abstract type AbstractImpurityModel <: AbstractImpurityHamiltonian end
sys_size(x::AbstractImpurityModel) = error("sys_size not implemented for model type $(typeof(x))")

hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityModel) = error("hybriddynamics not implemented for model $(typeof(model))")
hybriddynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; kwargs...) = hybriddynamics(vacuumstate(lattice), lattice, model; kwargs...)

sysdynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = sysdynamics!(copy(gmps), lattice, model; kwargs...)
sysdynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamics! not implemented for model $(typeof(model))")
sysdynamicsstepper!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamicsstepper! not implemented for model $(typeof(model))")
sysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = sysdynamics!(vacuumstate(lattice), lattice, model; kwargs...)


function sysdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, model::AbstractImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
	return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
end 


function sysdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward!(gmps, lattice, model; trunc=trunc) : sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	end
end 

function sysdynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice, model::AbstractImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
	else
		if branch == :+
			return sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		elseif branch == :-
			return sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
		end
	end
end 



sysdynamics_forward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamics_forward! not implemented for model $(typeof(model))")
sysdynamics_backward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamics_backward! not implemented for model $(typeof(model))")
sysdynamics_imaginary!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamics_imaginary! not implemented for model $(typeof(model))")


systhermalstate(lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = systhermalstate!(vacuumstate(lattice), lattice, model; kwargs...)
"""
	systhermalstate!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...)
	Set the initial state of the impurity to be a locally thermal state
"""
function systhermalstate!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...)
	return systhermalstate_iterative!(gmps, lattice, model; kwargs...)
end 

include("boundarycondition.jl")
include("accsysdynamics/accsysdynamics.jl")
include("sysinitstate.jl")

# predefined models
include("predefined/siam.jl")
include("predefined/irlm.jl")
include("predefined/skmodel.jl")

# general model
include("generalimpurity/generalimpurity.jl")