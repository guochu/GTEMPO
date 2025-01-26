# interface
# only the impurity Hamiltonian, no bath
abstract type AbstractImpurityHamiltonian end
# abstract type AbstractImpurityModel <: AbstractImpurityHamiltonian end
# sys_size(x::AbstractImpurityModel) = error("sys_size not implemented for model type $(typeof(x))")

# hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityModel) = error("hybriddynamics not implemented for model $(typeof(model))")
# hybriddynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; kwargs...) = hybriddynamics(vacuumstate(lattice), lattice, model; kwargs...)

"""
	sysdynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...)

Return the bare impurity dynamics as a GMPS, it is assumed that 
the build connection terms ⟨āᵢaⱼ⟩ (j=i-1 for 1 ≤ i ≤ N) is absorbed 
into the bare impurity dynamics, see the explicit treatment of the 
single-orbital Anderson model (SIAM) as an example

However, the boundary connection term with i=N, j=1, is independently 
treated in the function boundarycondition, since sometimes we have 
anti-periodic condition, and one may not want to directly perform this 
term, but calculate it as the summation of two GMPSs on the fly [see PRB 109, 165113 (2024)]

The initial state of the impurity is vacuum, if one would like to 
use a thermal state, simply apply systhermalstate on top of this function
"""
sysdynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = sysdynamics!(copy(gmps), lattice, model; kwargs...)
sysdynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = sysdynamics!(vacuumstate(lattice), lattice, model; kwargs...)
# sysdynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamics! not implemented for model $(typeof(model))")
sysdynamicsstepper!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamicsstepper! not implemented for model $(typeof(model))")

"""
	sysdynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...)

The inplace version of the function sysdynamics, which applys the bare impurity onto the given GMPS
"""
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

Reset the initial state of the impurity to be a local thermal state
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