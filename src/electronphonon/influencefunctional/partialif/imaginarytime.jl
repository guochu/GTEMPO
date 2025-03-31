

# naive implementation with N^2 gate operations
function hybriddynamics_naive!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	if lattice.bands == 1
		return hybriddynamics_1band_naive!(gmps, lattice, corr, trunc=trunc)
	else
		return hybriddynamics_2band_naive!(gmps, lattice, corr, trunc=trunc)
	end
end

function hybriddynamics_1band_naive!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr::ImagCorrelationFunction; kwargs...)
	@assert lattice.bands == 1
	return _hybriddynamics_1band_naive!(gmps, lattice, corr, 1; kwargs...)
end 

function _hybriddynamics_1band_naive!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr1::ImagCorrelationFunction, band::Int; trunc::TruncationScheme=DefaultITruncation)
	corr = corr1.data
	k = lattice.N
	for i in 1:k, j in 1:k
		pos1, pos2 = index(lattice, i, band=band), index(lattice, j, band=band)
		# coef = corr[i, j]
		coef = exp(corr[i, j]) - 1
		if i == j
			t = exp(NTerm(pos1, coeff=coef))
			apply!(t, gmps)
		else
			t = exp(NTerm(pos1, pos2, coeff=coef))
			apply!(t, gmps)
		end
		
		canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
		# println(bond_dimension(gmps), " ", norm(gmps), " ", coef)
	end
	return gmps
end

function hybriddynamics_2band_naive!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	@assert lattice.bands == 2
	for band in 1:lattice.bands
		_hybriddynamics_1band_naive!(gmps, lattice, corr1, band, trunc=trunc)
	end

	corr = corr1.data
	k = lattice.N
	for i in 1:k, j in 1:k
		pos1, pos2 = index(lattice, i, band=1), index(lattice, j, band=2)
		# coef = corr[i, j] + corr[j, i]
		# coef = exp(corr[i, j]) + exp(corr[j, i]) - 2
		coef = exp(corr[i, j] + corr[j, i]) - 1
		t = exp(NTerm(pos1, pos2, coeff=coef))
		apply!(t, gmps)
		canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
	end
	return gmps	
end


