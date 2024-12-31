"""
	hybriddynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, corr::ImagCorrelationFunction; band, trunc)

imaginary-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	corr = corr1.data
	k = lattice.k-1
	for i in 1:k
		tmp = partialif_hybrid(lattice, i+1, [0; view(corr, i, 1:k)], band=band)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end


function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	corr = corr1.data
	k = lattice.k-1
	for i in 1:k
		tmp = partialif_hybrid_naive(lattice, i+1, [0; view(corr, i, 1:k)], band=band, trunc=trunc)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end

function partialif_hybrid(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1)
	row = index(lattice, i, band=band, conj=true)
	col_pos = [index(lattice, j, band=band, conj=false) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end

function partialif_hybrid_naive(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	gmps = vacuumstate(lattice)
	for j in 1:lattice.k
		pos1, pos2 = index(lattice, i, conj=true, band=band), index(lattice, j, conj=false, band=band)
		t = exp(GTerm(pos1, pos2, coeff=cols[j]))
		apply!(t, gmps)
		canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
	end
	return gmps
end