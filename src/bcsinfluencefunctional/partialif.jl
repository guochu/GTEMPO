# column legs of a BCS partial IF: all branches and both flavor indices of the
# column operator, with `band2′` selected by the column conjugation flag c2
function _bcscolumns(lattice::AbstractGrassmannLattice, band1::Int, band2::Int)
	cols = GrassmannContourIndex[]
	for b2 in branches(lattice), c2 in (true, false)
		k2 = (b2 === :τ) ? lattice.Nτ : lattice.Nt
		band2′ = ifelse(c2, band2, band1)
		for j in 1:k2
			push!(cols, GrassmannContourIndex(j, b2, band2′, c2))
		end
	end
	return cols
end


# ---------------------------- contour-independent construction ----------------------------
"""
	hybriddynamics!(gmps, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; orbital=1, trunc=Defaults.alg_trunc_if())

Build the BCS partial IF on any (imaginary-, real- or mixed-time) Grassmann
lattice in place on `gmps`. For every branch the rows run over the
correlation-function nodes (`kτ - 1` for `:τ`, `kt - 1` for the real branches).
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; 
									orbital::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= orbital <= lattice.orbitals) || throw(BoundsError(1:lattice.orbitals, orbital))
	band1, band2 = 2*orbital-1, 2*orbital
	cols = _bcscolumns(lattice, band1, band2)
	for b1 in branches(lattice), c1 in (true, false)
		k1 = (b1 === :τ) ? lattice.kτ - 1 : lattice.kt - 1
		band1′ = ifelse(c1, band1, band2)
		for i in 1:k1
			tmp = partialif_hybrid(lattice, GrassmannContourIndex(i, b1, band1′, c1), cols, corr)
			gmps = mult!(gmps, tmp, trunc=trunc)[1]
		end
	end
	return gmps
end


hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction, alg::PartialIF; orbital::Int=1) = hybriddynamics(
				gmps, lattice, corr; orbital=orbital, trunc=alg.trunc)
hybriddynamics(lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction, alg::PartialIF; orbital::Int=1) = hybriddynamics!(
				vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; orbital=orbital, trunc=alg.trunc)

hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics(lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; kwargs...) = hybriddynamics!(vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; kwargs...)

hybriddynamics_naive(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; kwargs...) = hybriddynamics_naive!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics_naive(lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; kwargs...) = hybriddynamics_naive!(vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice), lattice, corr; kwargs...)


function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; 
									orbital::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= orbital <= lattice.orbitals) || throw(BoundsError(1:lattice.orbitals, orbital))
	band1, band2 = 2*orbital-1, 2*orbital
	cols = _bcscolumns(lattice, band1, band2)
	for b1 in branches(lattice), c1 in (true, false)
		k1 = (b1 === :τ) ? lattice.Nτ : lattice.Nt
		band1′ = ifelse(c1, band1, band2)
		for i in 1:k1
			tmp = partialif_hybrid_naive(lattice, GrassmannContourIndex(i, b1, band1′, c1), cols, corr, trunc=trunc)
			gmps = mult!(gmps, tmp, trunc=trunc)[1]
		end
	end
	return gmps
end
