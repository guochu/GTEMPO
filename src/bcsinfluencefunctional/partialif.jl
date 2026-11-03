# column legs of a BCS partial IF: for every branch and time, the three column
# field types (a↑, a↓, ā↓) covered by the exponent's second slots
# 列腿覆盖指数 (ā↑,a↑)/(ā↓,a↓)/(ā↑,ā↓)/(a↓,a↑) 的全部第二字段：
# (band1,false)=a↑、(band2,false)=a↓、(band2,true)=ā↓（见 partialif_hybrid 的
# BCS 系数分派；行腿自身的同节点列在装配处剔除）
function _bcscolumns(lattice::AbstractGrassmannLattice, band1::Int, band2::Int)
	cols = GrassmannContourIndex[]
	for b2 in branches(lattice)
		# 与 normal 的 `_columns` 一致：τ 列到 Nτ，实时列到 kt（全部节点）
		k2 = (b2 === :τ) ? lattice.Nτ : lattice.kt
		for j in 1:k2
			push!(cols, GrassmannContourIndex(j, b2, band1, false))
			push!(cols, GrassmannContourIndex(j, b2, band2, false))
			push!(cols, GrassmannContourIndex(j, b2, band2, true))
		end
	end
	return cols
end


# ---------------------------- contour-independent construction ----------------------------
"""
	hybriddynamics!(gmps, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; orbital=1, trunc=Defaults.alg_trunc_if())

Build the BCS partial IF on any (imaginary-, real- or mixed-time) Grassmann
lattice in place on `gmps`. For every branch the rows run over the
correlation-function nodes (`kτ - 1` for `:τ`, `kt` for the real branches),
matching the normal-bath construction.
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction;
									orbital::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= orbital <= lattice.orbitals) || throw(BoundsError(1:lattice.orbitals, orbital))
	band1, band2 = 2*orbital-1, 2*orbital
	cols = _bcscolumns(lattice, band1, band2)
	# 行腿 = 指数各项的第一字段（见 partialif_hybrid 的 BCS 系数分派）：
	# ā↑（正常 ↑↑ + 反常 ↑↓）、ā↓（正常 ↓↓）、a↓（反常 ↓↑）
	for b1 in branches(lattice), (band1′, c1) in ((band1, true), (band2, true), (band2, false))
		k1 = (b1 === :τ) ? lattice.kτ - 1 : lattice.kt
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
	for b1 in branches(lattice), (band1′, c1) in ((band1, true), (band2, true), (band2, false))
		k1 = (b1 === :τ) ? lattice.Nτ : lattice.kt
		for i in 1:k1
			tmp = partialif_hybrid_naive(lattice, GrassmannContourIndex(i, b1, band1′, c1), cols, corr, trunc=trunc)
			gmps = mult!(gmps, tmp, trunc=trunc)[1]
		end
	end
	return gmps
end
