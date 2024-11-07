
# function retardedinteractdynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::AbstractMixedCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
# 	if LayoutStyle(lattice) isa BandLocalLayout
# 		return _retardedinteractdynamics!(gmps, lattice, corr, trunc=trunc)
# 	else
# 		lattice2 = similar(lattice, ordering=A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1())
# 		gmps2 = _retardedinteractdynamics!(vacuumstate(lattice2), lattice2, corr; trunc=trunc)
# 		perm = matchindices2(lattice, lattice2)
# 		gmps2 = permute(gmps2, perm, trunc=trunc)
# 		return mult!(gmps, gmps2, trunc=trunc)
# 	end
# end


function retardedinteractdynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::AbstractMixedCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	if lattice.bands == 1
		return retardedinteractdynamics_1band!(gmps, lattice, corr, trunc=trunc)
	else
		return retardedinteractdynamics_2band!(gmps, lattice, corr, trunc=trunc)
	end
end


function retardedinteractdynamics_1band!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::AbstractMixedCorrelationFunction; 
											trunc::TruncationScheme=DefaultITruncation)	
	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	@assert lattice.bands == 1
	band = 1
	for i in 1:lattice.kt-1, b1 in (:+, :-)
		pos1a, pos1b = index(lattice, i+1, band=band, conj=true, branch=b1), index(lattice, i, band=band, conj=false, branch=b1)
		for j in 1:lattice.kt-1, b2 in (:+, :-)
			pos2a, pos2b = index(lattice, j+1, band=band, conj=true, branch=b2), index(lattice, j, band=band, conj=false, branch=b2)
			if (i == j) && (b1 == b2)
				t = exp(GTerm(pos1a, pos1b, coeff=index(corr, i, j, b1=b1, b2=b2)))
			else
				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=index(corr, i, j, b1=b1, b2=b2)))
			end
			apply!(t, gmps)
			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
		end
		for j in 1:lattice.kτ-1
			b2 = :τ
			pos2a, pos2b = index(lattice, j+1, band=band, conj=true, branch=b2), index(lattice, j, band=band, conj=false, branch=b2)
			t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=index(corr, i, j, b1=b1, b2=b2)))
			apply!(t, gmps)
			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
		end
	end
	for i in 1:lattice.kτ-1
		b1 = :τ
		pos1a, pos1b = index(lattice, i+1, band=band, conj=true, branch=b1), index(lattice, i, band=band, conj=false, branch=b1)
		for j in 1:lattice.kt-1, b2 in (:+, :-)
			pos2a, pos2b = index(lattice, j+1, band=band, conj=true, branch=b2), index(lattice, j, band=band, conj=false, branch=b2)
			t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=index(corr, i, j, b1=b1, b2=b2)))
			apply!(t, gmps)
			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
		end
		for j in 1:lattice.kτ-1
			b2 = :τ
			pos2a, pos2b = index(lattice, j+1, band=band, conj=true, branch=b2), index(lattice, j, band=band, conj=false, branch=b2)
			if (i == j)
				t = exp(GTerm(pos1a, pos1b, coeff=index(corr, i, j, b1=b1, b2=b2)))
			else
				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=index(corr, i, j, b1=b1, b2=b2)))
			end
			apply!(t, gmps)
			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
		end
	end
	return gmps
end


function retardedinteractdynamics_2band!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::AbstractMixedCorrelationFunction; 
											trunc::TruncationScheme=DefaultMPOTruncation)
	@assert lattice.bands == 2
	error("not implemented")
end

