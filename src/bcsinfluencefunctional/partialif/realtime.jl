

# function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::BCSCorrelationFunction{<:RealCorrelationFunction}; 
# 									orbital::Int=1, trunc::TruncationScheme=DefaultITruncation)
# 	(1 <= orbital <= lattice.orbitals) || throw(BoundsError(1:lattice.orbitals, orbital))
# 	orth = Orthogonalize(TK.SVD(), trunc)
# 	band1, band2 = 2*orbital-1, 2*orbital
# 	k = lattice.k-1
# 	for i in 1:k, b1 in branches(lattice), c1 in (true, false)
# 		tmp = vacuumstate(lattice)
# 		band1′ = ifelse(c1, band1, band2)
# 		pos1 = index(lattice, i, conj=c1, band=band1′, branch=b1)
# 		for j in 1:k, b2 in branches(lattice), c2 in (true, false)
# 			coef = index(corr[c1, c2], i, j, b1=b1, b2=b2)
# 			band2′ = ifelse(c2, band2, band1)
# 			pos2 = index(lattice, j, conj=c2, band=band2′, branch=b2)
# 			t = exp(GTerm(pos1, pos2, coeff=coef))
# 			apply!(t, tmp)
# 			canonicalize!(tmp, alg=orth)
# 		end
# 		gmps = mult!(gmps, tmp, trunc=trunc)
# 	end
# 	return gmps
# end


function hybriddynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::BCSCorrelationFunction{<:RealCorrelationFunction}; 
									orbital::Int=1, trunc::TruncationScheme=DefaultITruncation)
	(1 <= orbital <= lattice.orbitals) || throw(BoundsError(1:lattice.orbitals, orbital))
	orth = Orthogonalize(TK.SVD(), trunc)
	band1, band2 = 2*orbital-1, 2*orbital
	k = lattice.k-1
	for i in 1:k, b1 in branches(lattice), c1 in (true, false)
		band1′ = ifelse(c1, band1, band2)
		pos1 = index(lattice, i, conj=c1, band=band1′, branch=b1)
		pos2s = Int[]
		coefs = scalartype(lattice)[]
		for j in 1:k, b2 in branches(lattice), c2 in (true, false)
			coef = index(corr[c1, c2], i, j, b1=b1, b2=b2)
			band2′ = ifelse(c2, band2, band1)
			pos2 = index(lattice, j, conj=c2, band=band2′, branch=b2)
			push!(pos2s, pos2)
			push!(coefs, coef)
		end
		tmp = partialmpo(pos1, pos2s, coefs) * vacuumstate(lattice)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end