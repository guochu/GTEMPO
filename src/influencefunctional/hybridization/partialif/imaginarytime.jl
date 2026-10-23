"""
	hybriddynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, corr::ImagCorrelationFunction; band, trunc)

imaginary-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	k = lattice.k-1
	for i in 1:k
		row = [corr1[GrassmannContourIndex(i, :τ, band, false), GrassmannContourIndex(j, :τ, band, true)] for j in 1:k]
		tmp = partialif_hybrid(lattice, i+1, [0; row], band=band)
		gmps = mult!(gmps, tmp, trunc=trunc)[1]
	end
	return gmps
end


function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	k = lattice.k-1
	for i in 1:k
		row = [corr1[GrassmannContourIndex(i, :τ, band, false), GrassmannContourIndex(j, :τ, band, true)] for j in 1:k]
		tmp = partialif_hybrid_naive(lattice, i+1, [0; row], band=band, trunc=trunc)
		gmps = mult!(gmps, tmp, trunc=trunc)[1]
	end
	return gmps
end

"""
	partialif_hybrid(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band)

Build the partial IF I = e^{āᵢΣⱼΔᵢⱼaⱼ} analytically as an MPO with bond dimension 2
cols: the i-th row of the hybridization function Δᵢⱼ
"""
function partialif_hybrid(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1)
	row = lattice[GrassmannContourIndex(i, :τ, band, true)]
	col_pos = [lattice[GrassmannContourIndex(j, :τ, band, false)] for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(promote_type(scalartype(lattice), eltype(cols)), lattice)
end

function partialif_hybrid_naive(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	gmps = vacuumstate(promote_type(scalartype(lattice), eltype(cols)), lattice)
	for j in 1:lattice.k
		pos1, pos2 = lattice[GrassmannContourIndex(i, :τ, band, true)], lattice[GrassmannContourIndex(j, :τ, band, false)]
		t = exp(GTerm(pos1, pos2, coeff=cols[j]))
		apply!(t, gmps)
		canonicalize!(gmps, alg=Orthogonalize(orth=TK.SVD(), trunc=trunc))
	end
	return gmps
end