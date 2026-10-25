include("util.jl")

# static mode, compute the whole IF and then measure

"""
	hybriddynamics!(gmps::FockMPS, lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; trunc)

Build the phonon (retarded-interaction) influence functional
`exp(Σᵢⱼ Δᵢⱼ N̂ᵢ N̂ⱼ)` on a Fock lattice in place on `gmps`, where
`N̂ᵢ = Σ_band n̂ᵢ` is the total occupation of contour node `i`.

The IF is the product of one partial IF per contour leg `(time, branch, band)`,
each built by [`partialif_hybrid`](@ref); the contour branches (`:τ`, `:+`,
`:-`) and the bands are read from `lattice`, so the same implementation covers
imaginary, real and mixed time (no per-lattice methods).
"""
function hybriddynamics!(gmps::FockMPS, lattice::AbstractFockLattice, corr::AbstractCorrelationFunction;
						 trunc::TruncationScheme=Defaults.alg_trunc_if())
	cols = _fockcolumns(lattice)
	for row in cols
		mult!(gmps, partialif_hybrid(lattice, row, cols, corr), trunc=trunc)
	end
	return gmps
end

"""
	hybriddynamics_naive!(gmps::FockMPS, lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; trunc)

Naive `N²`-gate construction of the same phonon IF: apply each
`exp(Δᵢⱼ n̂ᵢ n̂ⱼ)` as a two-body `ExpNTerm` (one-body for `i == j`) and
canonicalize after every step. Works for any number of bands and all three
contours.
"""
function hybriddynamics_naive!(gmps::FockMPS, lattice::AbstractFockLattice, corr::AbstractCorrelationFunction;
							   trunc::TruncationScheme=Defaults.alg_trunc_if())
	cols = _fockcolumns(lattice)
	alg = Orthogonalize(orth=TK.SVD(), trunc=trunc)
	for row in cols
		pos1 = index(lattice, row.j, branch=row.branch, band=row.band)
		for col in cols
			pos2 = index(lattice, col.j, branch=col.branch, band=col.band)
			# ExpNTerm 的系数是 Λ = e^{Δ}（占位 n̂ᵢn̂ⱼ = 1 时的因子）
			c = exp(corr[row, col])
			t = pos1 == pos2 ? ExpNTerm(pos1, coeff=c) : ExpNTerm(pos1, pos2, coeff=c)
			apply!(t, gmps)
			canonicalize!(gmps, alg=alg)
		end
	end
	return gmps
end

hybriddynamics(gmps::FockMPS, lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics(lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(
				vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; kwargs...)


hybriddynamics_naive(gmps::FockMPS, lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics_naive(lattice::AbstractFockLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics_naive!(
				vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; kwargs...)
