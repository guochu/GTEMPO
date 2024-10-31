hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::PartialIF; band::Int=1) = hybriddynamics(
				gmps, lattice, corr; band=band, trunc=alg.trunc)
hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::PartialIF; band::Int=1) = hybriddynamics!(
				vacuumstate(lattice), lattice, corr; band=band, trunc=alg.trunc)

hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(vacuumstate(lattice), lattice, corr; kwargs...)


### for single impurity models

"""
	qim_hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::NTuple{4, <:AbstractMatrix}; trunc)

Real time hybrid dynamics for single impurity model, all the bands 
share the same bath 
"""
function qim_hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...)
	for band in 1:lattice.bands
		gmps = hybriddynamics!(gmps, lattice, corr; band=band, kwargs...)
	end
	return gmps
end
qim_hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
qim_hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamics!(vacuumstate(lattice), lattice, corr; kwargs...)

