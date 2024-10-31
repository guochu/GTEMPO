hybriddynamicsstepper(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamicsstepper!(copy(gmps), lattice, corr; kwargs...)


# function qim_hybriddynamicsstepper!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...)
# 	for band in 1:lattice.bands
# 		gmps = hybriddynamicsstepper!(gmps, lattice, corr; band=band, kwargs...)
# 	end
# 	return gmps
# end
# qim_hybriddynamicsstepper(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamicsstepper!(copy(gmps), lattice, corr; kwargs...)
