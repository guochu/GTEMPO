include("util.jl")

# ---------------------------------------------------------------------------
# static mode: compute the whole IF and then measure
# ---------------------------------------------------------------------------

"""
	hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg; band::Int=1) 

Build the IF=e^{ΣᵢⱼΔᵢⱼāᵢaⱼ} as a GMPS using algorithm alg
corr: fermionic hybridization function calculated using QuAPI
"""
hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::PartialIF; band::Int=1) = hybriddynamics(
				gmps, lattice, corr; band=band, trunc=alg.trunc)
hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::PartialIF; band::Int=1) = hybriddynamics!(
				vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; band=band, trunc=alg.trunc)

hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(
				vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; kwargs...)

hybriddynamics_naive(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics_naive(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics_naive!(
				vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; kwargs...)

"""
	hybriddynamics!(gmps, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; band=1, trunc=Defaults.alg_trunc_if())
	hybriddynamics_naive!(gmps, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; band=1, trunc=Defaults.alg_trunc_if())

Build the MPS-IF of the hybridization function in place on `gmps`, for a
single band. The rows are the contour legs `āᵢ` of every branch of `lattice`
(imaginary-, real- or mixed-time), and each row is applied as one partial IF
built by [`partialif_hybrid`](@ref) (or [`partialif_hybrid_naive`](@ref)).
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	cols = _columns(lattice, band)
	for b1 in branches(lattice)
		k1 = (b1 === :τ) ? lattice.Nτ : lattice.kt
		for i in 1:k1
			tmp = partialif_hybrid(lattice, GrassmannContourIndex(i, b1, band, true), cols, corr)
			gmps = mult!(gmps, tmp, trunc=trunc)[1]
		end
	end
	return gmps
end

function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	cols = _columns(lattice, band)
	for b1 in branches(lattice)
		k1 = (b1 === :τ) ? lattice.Nτ : lattice.kt
		for i in 1:k1
			tmp = partialif_hybrid_naive(lattice, GrassmannContourIndex(i, b1, band, true), cols, corr, trunc=trunc)
			gmps = mult!(gmps, tmp, trunc=trunc)[1]
		end
	end
	return gmps
end


# ### for single impurity models

# """
# 	qim_hybriddynamics(gmps::GrassmannMPS, lattice::AbstractRealGrassmannLattice, corr::NTuple{4, <:AbstractMatrix}; trunc)

# Real time hybrid dynamics for single impurity model, all the bands 
# share the same bath 
# """
# function qim_hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...)
# 	for band in 1:lattice.bands
# 		gmps = hybriddynamics!(gmps, lattice, corr; band=band, kwargs...)
# 	end
# 	return gmps
# end
# qim_hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
# qim_hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamics!(vacuumstate(lattice), lattice, corr; kwargs...)


# ---------------------------------------------------------------------------
# dynamic mode: evolve the IF and measure during the evolution, similar to TEMPO
# ---------------------------------------------------------------------------
include("realtime_stepper.jl")
