function hybriddynamics!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	for i in 1:lattice.Nt, b1 in branches(lattice), band1 in 1:lattice.bands
		pos1 = index(lattice, i, band=band1, branch=b1)
		pos2s = Int[]
		coefs = scalartype(lattice)[]
		for j in 1:lattice.Nt, b2 in branches(lattice), band2 in 1:lattice.bands
			pos2 = index(lattice, j, band=band2, branch=b2)
			coef = index(corr, i, j, b1=b1, b2=b2)
			push!(pos2s, pos2)
			push!(coefs, coef)
		end
		tmp = partialif_densemps(length(lattice), pos1, pos2s, coefs)
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end


# function hybriddynamics!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
# 	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
# 	if lattice.bands == 1
# 		return hybriddynamics_1band!(gmps, lattice, corr, trunc=trunc)
# 	else
# 		return hybriddynamics_2band!(gmps, lattice, corr, trunc=trunc)
# 	end
# end

# function hybriddynamics_1band!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction; kwargs...)
# 	@assert lattice.bands == 1
# 	return _hybriddynamics_1band!(gmps, lattice, corr, 1; kwargs...)
# end 


# function _hybriddynamics_1band!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction, band::Int; 
# 												trunc::TruncationScheme=DefaultITruncation)	
# 	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
# 	alg = Orthogonalize(TK.SVD(), trunc)
# 	for i in 1:lattice.Nt, b1 in branches(lattice)
# 		pos1 = index(lattice, i, band=band, branch=b1)
# 		# tmp = vacuumstate(lattice)
# 		# for j in 1:lattice.Nt, b2 in branches(lattice)
# 		# 	pos2 = index(lattice, j, band=band, branch=b2)
# 		# 	c = exp(index(corr, i, j, b1=b1, b2=b2)) 
# 		# 	if pos1 == pos2
# 		# 		t = ExpNTerm(pos1, coeff=c)
# 		# 	else
# 		# 		t = ExpNTerm(pos1, pos2, coeff=c)
# 		# 	end
# 		# 	apply!(t, tmp)
# 		# 	canonicalize!(tmp, alg=alg)
# 		# end
# 		pos2s = Int[]
# 		coefs = scalartype(lattice)[]
# 		for j in 1:lattice.Nt, b2 in branches(lattice)
# 			pos2 = index(lattice, j, band=band, branch=b2)
# 			c = index(corr, i, j, b1=b1, b2=b2)
# 			push!(pos2s, pos2)
# 			push!(coefs, c)
# 		end
# 		tmp = partialif_densemps(length(lattice), pos1, pos2s, coefs)
# 		mult!(gmps, tmp, trunc=trunc)
# 	end
# 	return gmps
# end

# function hybriddynamics_2band!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultMPOTruncation)
# 	@assert lattice.bands == 2
# 	_hybriddynamics_1band!(gmps, lattice, corr, 2, trunc=trunc)
# 	alg = Orthogonalize(TK.SVD(), trunc)
# 	for i in 1:lattice.Nt, b1 in branches(lattice)
# 		pos1 = index(lattice, i, band=1, branch=b1)
# 		tmp = vacuumstate(lattice)
# 		# for j in 1:lattice.Nt, b2 in branches(lattice), band2 in 1:lattice.bands
# 		# 	pos2 = index(lattice, j, band=band2, branch=b2)
# 		# 	c = index(corr, i, j, b1=b1, b2=b2)
# 		# 	coef = ifelse(band2==1, exp(c) , exp(c + index(corr, j, i, b1=b2, b2=b1)) )
# 		# 	if pos1 == pos2
# 		# 		t = ExpNTerm(pos1, coeff=coef)
# 		# 	else
# 		# 		t = ExpNTerm(pos1, pos2, coeff=coef)
# 		# 	end
# 		# 	apply!(t, tmp)
# 		# 	canonicalize!(tmp, alg=alg)
# 		# end
# 		pos2s = Int[]
# 		coefs = scalartype(lattice)[]
# 		for j in 1:lattice.Nt, b2 in branches(lattice), band2 in 1:lattice.bands
# 			pos2 = index(lattice, j, band=band2, branch=b2)
# 			c = index(corr, i, j, b1=b1, b2=b2)
# 			coef = ifelse(band2==1, c , c + index(corr, j, i, b1=b2, b2=b1))
# 			push!(pos2s, pos2)
# 			push!(coefs, coef)
# 		end
# 		tmp = partialif_densemps(length(lattice), pos1, pos2s, coefs)
# 		mult!(gmps, tmp, trunc=trunc)
# 	end
# 	return gmps
# end


# naive implementation with N^2 gate operations
function hybriddynamics_naive!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	if lattice.bands == 1
		return hybriddynamics_1band_naive!(gmps, lattice, corr, trunc=trunc)
	else
		return hybriddynamics_2band_naive!(gmps, lattice, corr, trunc=trunc)
	end
end

function hybriddynamics_1band_naive!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction; kwargs...)
	@assert lattice.bands == 1
	return _hybriddynamics_1band_naive!(gmps, lattice, corr, 1; kwargs...)
end 


function _hybriddynamics_1band_naive!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction, band::Int; 
												trunc::TruncationScheme=DefaultITruncation)	
	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	alg = Orthogonalize(TK.SVD(), trunc)
	for b1 in branches(lattice), b2 in branches(lattice)
		for i in 1:lattice.Nt, j in 1:lattice.Nt
			pos1, pos2 = index(lattice, i, band=band, branch=b1), index(lattice, j, band=band, branch=b2)
			c = exp(index(corr, i, j, b1=b1, b2=b2)) 
			if pos1 == pos2
				t = ExpNTerm(pos1, coeff=c)
			else
				t = ExpNTerm(pos1, pos2, coeff=c)
			end
			apply!(t, gmps)
			canonicalize!(gmps, alg=alg)
		end
	end
	return gmps
end

function hybriddynamics_2band_naive!(gmps::FockMPS, lattice::RealFockLattice, corr::RealCorrelationFunction; 
											trunc::TruncationScheme=DefaultMPOTruncation)
	@assert lattice.bands == 2
	for band in 1:lattice.bands
		_hybriddynamics_1band_naive!(gmps, lattice, corr, band, trunc=trunc)
	end
	alg = Orthogonalize(TK.SVD(), trunc)
	for i in 1:lattice.Nt, b1 in branches(lattice)
		pos1 = index(lattice, i, band=1, branch=b1)
		for j in 1:lattice.Nt, b2 in branches(lattice)
			pos2 = index(lattice, j, band=2, branch=b2)
			c = exp(index(corr, i, j, b1=b1, b2=b2) + index(corr, j, i, b1=b2, b2=b1)) 
			t = ExpNTerm(pos1, pos2, coeff=c)
			apply!(t, gmps)
			canonicalize!(gmps, alg=alg)
		end
	end
	return gmps
end


