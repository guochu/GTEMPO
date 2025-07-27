include("partialif/partialif.jl")




function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::BCSCorrelationFunction; 
									orbital::Int=1, trunc::TruncationScheme=DefaultITruncation)
	(1 <= orbital <= lattice.orbitals) || throw(BoundsError(1:lattice.orbitals, orbital))
	orth = Orthogonalize(TK.SVD(), trunc)
	band1, band2 = 2*orbital-1, 2*orbital
	for b1 in branches(lattice), c1 in (true, false)
		k1 = (b1 == :τ) ? lattice.Nτ : lattice.Nt
		for i in 1:k1
			tmp = vacuumstate(scalartype(corr), lattice)
			band1′ = ifelse(c1, band1, band2)
			i′ = (b1 == :τ) ? i+1 : i
			pos1 = index(lattice, i′, conj=c1, band=band1′, branch=b1)
			for b2 in branches(lattice), c2 in (true, false)
				k2 = (b2 == :τ) ? lattice.Nτ : lattice.Nt
				# _corr = branch(corr[c1, c2], b1=b1, b2=b2)
				for j in 1:k2
					coef = index(corr[c1, c2], i, j, b1=b1, b2=b2)
					# coef = _corr[i, j]
					band2′ = ifelse(c2, band2, band1)
					j′ = (b1==:τ) ? j+1 : j
					pos2 = index(lattice, j′, conj=c2, band=band2′, branch=b2)
					t = exp(GTerm(pos1, pos2, coeff=coef))
					apply!(t, tmp)
					canonicalize!(tmp, alg=orth)
				end
			end
			gmps = mult!(gmps, tmp, trunc=trunc)			
		end
	end
	return gmps
end