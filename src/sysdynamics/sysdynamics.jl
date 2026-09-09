# interface
# only the impurity Hamiltonian, no bath
abstract type AbstractImpurityHamiltonian end
# constant Hamiltonian (quench hamiltonian belongs to ConstImpurityHamiltonian)
# Only ConstImpurityHamiltonian supports sysdynamics_fast 
abstract type ConstImpurityHamiltonian <: AbstractImpurityHamiltonian end
# time-dependent Hamiltonian
abstract type GeneralTdImpurityHamiltonian <: AbstractImpurityHamiltonian end

# interface for AbstractImpurityHamiltonian 
fock_propagator(model::ConstImpurityHamiltonian, branch::Symbol, dt::Real) = error("fock_propagator not implemented for model $(typeof(model))")
fock_thermalstate(h::AbstractImpurityHamiltonian, β::Real) = error("fock_thermalstate not implemented for model $(typeof(h))")
fock_propagator(model::GeneralTdImpurityHamiltonian, branch::Symbol, dt::Real, t::Real) = error("fock_propagator not implemented for model $(typeof(model))")

num_bands(h::AbstractImpurityHamiltonian) = h.bands

# hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityModel) = error("hybriddynamics not implemented for model $(typeof(model))")
# hybriddynamics(lattice::AbstractGrassmannLattice, model::AbstractImpurityModel; kwargs...) = hybriddynamics(vacuumstate(lattice), lattice, model; kwargs...)

sysdynamicsstepper!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::ConstImpurityHamiltonian; kwargs...) = error("sysdynamicsstepper! not implemented for model $(typeof(model))")

# Deprecated Trotterized dynamics: the exact propagator sysdynamics /
# baresysdynamics (built from the Fock-space propagator) supersedes
# these; kept for reference and comparison
"""
	sysdynamics_deprecated(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...)

Deprecated Trotterized impurity dynamics, superseded by the exact
`sysdynamics`. The bulk connection terms ⟨āᵢaⱼ⟩ (j=i-1 for 1 ≤ i ≤ N)
are absorbed into the propagator; the boundary connection term with
i=N, j=1 is treated separately in `boundarycondition`.
"""
sysdynamics_deprecated(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = sysdynamics_deprecated!(copy(gmps), lattice, model; kwargs...)
sysdynamics_deprecated(lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = sysdynamics_deprecated!(vacuumstate(lattice), lattice, model; kwargs...)

"""
	sysdynamics_deprecated!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...)

The inplace version of `sysdynamics_deprecated`.
"""
function sysdynamics_deprecated!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, model::AbstractImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
	return sysdynamics_imaginary_deprecated!(gmps, lattice, model; trunc=trunc)
end


function sysdynamics_deprecated!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian;
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward_deprecated!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_backward_deprecated!(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward_deprecated!(gmps, lattice, model; trunc=trunc) : sysdynamics_backward_deprecated!(gmps, lattice, model; trunc=trunc)
	end
end

function sysdynamics_deprecated!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice, model::AbstractImpurityHamiltonian;
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		sysdynamics_forward_deprecated!(gmps, lattice, model; trunc=trunc)
		sysdynamics_backward_deprecated!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_imaginary_deprecated!(gmps, lattice, model; trunc=trunc)
	else
		if branch == :+
			return sysdynamics_forward_deprecated!(gmps, lattice, model; trunc=trunc)
		elseif branch == :-
			return sysdynamics_backward_deprecated!(gmps, lattice, model; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return sysdynamics_imaginary_deprecated!(gmps, lattice, model; trunc=trunc)
		end
	end
end



# per-model Taylor-expanded steps of the deprecated dynamics
sysdynamics_forward_deprecated!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamics_forward_deprecated! not implemented for model $(typeof(model))")
sysdynamics_backward_deprecated!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamics_backward_deprecated! not implemented for model $(typeof(model))")
sysdynamics_imaginary_deprecated!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) = error("sysdynamics_imaginary_deprecated! not implemented for model $(typeof(model))")

# predefined models
include("predefined/siam.jl")
include("predefined/irlm.jl")
include("predefined/skmodel.jl")

# general model
include("generalimpurity/generalimpurity.jl")
