"""
	hybriddynamics(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::MixedCorrelationFunction; band::Int, trunc)

real-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::AbstractMixedCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	kt, Nτ = lattice.kt, lattice.Nτ
	cols = _mixedcolumns(kt, Nτ, band)
	for b1 in (:+, :-)
		for i in 1:kt
			tmp = partialif_hybrid(lattice, GrassmannContourIndex(i, b1, band, false), cols, corr)
			gmps = mult!(gmps, tmp, trunc=trunc)[1]
		end
	end
	for i in 1:Nτ
		tmp = partialif_hybrid(lattice, GrassmannContourIndex(i, :τ, band, false), cols, corr)
		gmps = mult!(gmps, tmp, trunc=trunc)[1]
	end
	return gmps
end

function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::AbstractMixedCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	kt, Nτ = lattice.kt, lattice.Nτ
	cols = _mixedcolumns(kt, Nτ, band)
	for b1 in (:+, :-)
		for i in 1:kt
			tmp = partialif_hybrid_naive(lattice, GrassmannContourIndex(i, b1, band, false), cols, corr, trunc=trunc)
			gmps = mult!(gmps, tmp, trunc=trunc)[1]
		end
	end
	for i in 1:Nτ
		tmp = partialif_hybrid_naive(lattice, GrassmannContourIndex(i, :τ, band, false), cols, corr, trunc=trunc)
		gmps = mult!(gmps, tmp, trunc=trunc)[1]
	end
	return gmps		
end

# column indices shared by all rows of a mixed-contour partial IF: the
# imaginary-time columns first, then the real-time columns (+ and - per step)
function _mixedcolumns(kt::Int, Nτ::Int, band::Int)
	cols = GrassmannContourIndex[]
	sizehint!(cols, Nτ + 2*kt)
	for j in 1:Nτ
		push!(cols, GrassmannContourIndex(j, :τ, band, true))
	end
	for j in 1:kt
		push!(cols, GrassmannContourIndex(j, :+, band, true))
		push!(cols, GrassmannContourIndex(j, :-, band, true))
	end
	return cols
end

function partialif_hybrid(lattice::MixedGrassmannLattice, i::GrassmannContourIndex,
						  cols::AbstractVector{<:GrassmannContourIndex}, corr)
	@assert length(cols) == lattice.Nτ + 2*lattice.kt
	row = lattice[_latticeindex(i, conj=true)]
	col_pos = [lattice[_latticeindex(c, conj=false)] for c in cols]
	coefs = [corr[i, c] for c in cols]
	mpo = partialmpo(row, col_pos, coefs)
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(coefs)), lattice)
end

function partialif_hybrid_naive(lattice::MixedGrassmannLattice, i::GrassmannContourIndex,
						  cols::AbstractVector{<:GrassmannContourIndex}, corr;
						  trunc::TruncationScheme=Defaults.alg_trunc_if())
	@assert length(cols) == lattice.Nτ + 2*lattice.kt
	alg = Orthogonalize(orth=TK.SVD(), trunc=trunc)
	gmps = vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	pos1 = lattice[_latticeindex(i, conj=true)]
	for c in cols
		pos2 = lattice[_latticeindex(c, conj=false)]
		t = exp(GTerm(pos1, pos2, coeff=corr[i, c]))
		apply!(t, gmps)
		canonicalize!(gmps, alg=alg)
	end
	return gmps
end
