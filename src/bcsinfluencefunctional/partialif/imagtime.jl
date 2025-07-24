

function bcs_hybriddynamics_naive!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr::BCSCorrelationFunction{<:ImagCorrelationFunction}; 
									orb::Int=1, trunc::TruncationScheme=DefaultITruncation)
	(1 <= orb <= lattice.norb) || throw(BoundsError(1:lattice.norb, orb))
	orth = Orthogonalize(TK.SVD(), trunc)
	# Δ11, Δ12, Δ21, Δ22 = corr[1, 1], corr[1,2], corr[2,1], corr[2,2]
	b1, b2 = 2*orb-1, 2*orb
	k = lattice.k-1
	for i in 1:k, c1 in (true, false)
		tmp = vacuumstate(lattice)
		pos1 = index(lattice, i, conj=c1, band=b1)
		for j in 1:k, c2 in (true, false)
			coef = corr[c1, c2][i, j]
			b2′ = ifelse(c2==c1, b2, b1)
			pos2 = index(lattice, j, conj=c2, band=b2′)
			t = exp(GTerm(pos1, pos2, coeff=coef))
			apply!(t, tmp)
			canonicalize!(tmp, alg=orth)
		end
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end

