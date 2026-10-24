# ---------------------------------------------------------------------------
# one time step of the IF dynamics
# ---------------------------------------------------------------------------

"""
	hybriddynamicsstepper!(gmps, lattice, corr; band=1, trunc=Defaults.alg_trunc_if())

In-place single time step of the influence functional dynamics: applies the
IF terms connecting the current time step to all previous steps onto `gmps`
and moves it forward by one step. The non-mutating version
`hybriddynamicsstepper` operates on a copy.
"""
function hybriddynamicsstepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	@assert length(lattice) == length(gmps) 
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	if lattice.k == 2
		for b1 in (:+, :-), b2 in (:+, :-)
			pos1 = lattice[GrassmannContourIndex(1, b1, band, true)]
			pos2 = lattice[GrassmannContourIndex(1, b2, band, false)]
			apply!(exp(GTerm(pos1, pos2, coeff=branch(corr; b1=b1, b2=b2)[1,1])), gmps)
		end
		canonicalize!(gmps, alg=Orthogonalize(orth=SVD(), trunc=trunc))
	end
	k = lattice.k
	cols = [GrassmannContourIndex(j, b2, band, false) for j in 1:k for b2 in (:+, :-)]
	for b1 in (:+, :-)
		tmp = partialif_hybrid(lattice, GrassmannContourIndex(k, b1, band, true), cols, corr)
		gmps = mult!(gmps, tmp, trunc=trunc)[1]
	end
	rows = [GrassmannContourIndex(i, b1, band, true) for i in 1:k-1 for b1 in (:+, :-)]
	for b2 in (:+, :-)
		tmp = partialif_hybrid(lattice, rows, GrassmannContourIndex(k, b2, band, false), corr)
		gmps = mult!(gmps, tmp, trunc=trunc)[1]
	end
	return gmps
end


# ---------------------------------------------------------------------------
# non-mutating wrapper
# ---------------------------------------------------------------------------
hybriddynamicsstepper(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamicsstepper!(copy(gmps), lattice, corr; kwargs...)


# function qim_hybriddynamicsstepper!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...)
# 	for band in 1:lattice.bands
# 		gmps = hybriddynamicsstepper!(gmps, lattice, corr; band=band, kwargs...)
# 	end
# 	return gmps
# end
# qim_hybriddynamicsstepper(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamicsstepper!(copy(gmps), lattice, corr; kwargs...)
