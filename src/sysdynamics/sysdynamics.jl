"""
	AbstractImpurityHamiltonian

Supertype of all impurity Hamiltonians. Interface: `fock_propagator`
(both subtypes) and `fock_thermalstate` (needed to set the thermal
initial state on real-time lattices); the number of impurity bands is
reported by `num_bands`.
"""
abstract type AbstractImpurityHamiltonian end
"""
	ConstImpurityHamiltonian <: AbstractImpurityHamiltonian

Constant (time-independent) impurity Hamiltonian; quench Hamiltonians
belong here. Only `ConstImpurityHamiltonian` supports `sysdynamics_fast`.
Interface: `fock_propagator(model, branch, dt)`.
"""
abstract type ConstImpurityHamiltonian <: AbstractImpurityHamiltonian end
"""
	AbstractTdImpurityHamiltonian <: AbstractImpurityHamiltonian

Time-dependent impurity Hamiltonian; only `sysdynamics` is supported.
Interface: `fock_propagator(model, branch, dt, t)` on the real-time
branches (`:+`/`:-`), with `t` the left endpoint of the physical time
interval of the step.
"""
abstract type AbstractTdImpurityHamiltonian <: AbstractImpurityHamiltonian end

# interface for AbstractImpurityHamiltonian 
fock_propagator(model::ConstImpurityHamiltonian, branch::Symbol, dt::Real) = error("fock_propagator not implemented for model $(typeof(model))")
fock_thermalstate(h::AbstractImpurityHamiltonian, β::Real) = error("fock_thermalstate not implemented for model $(typeof(h))")
fock_propagator(model::AbstractTdImpurityHamiltonian, branch::Symbol, dt::Real, t::Real) = error("fock_propagator not implemented for model $(typeof(model))")
# constant models: the propagator on the real-time branches is time-independent,
# the branch time t is ignored
fock_propagator(model::ConstImpurityHamiltonian, branch::Symbol, dt::Real, t::Real) = fock_propagator(model, branch, dt)

"""
	num_bands(h) -> Int

Number of impurity bands of the model. Defaults to the `bands` field of
the model; models without such a field define their own method (e.g.
`AndersonIM` reports 2 and `ToulouseIM` reports 1).
"""
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
