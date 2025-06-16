
function hybriddynamics!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	if lattice.bands == 1
		return hybriddynamics_1band!(gmps, lattice, corr, trunc=trunc)
	else
		return hybriddynamics_2band!(gmps, lattice, corr, trunc=trunc)
	end
end

function hybriddynamics_1band!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr::ImagCorrelationFunction; kwargs...)
	@assert lattice.bands == 1
	return _hybriddynamics_1band!(gmps, lattice, corr, 1; kwargs...)
end 

function _hybriddynamics_1band!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr1::ImagCorrelationFunction, band::Int; trunc::TruncationScheme=DefaultITruncation)
	corr = corr1.data
	k = lattice.N
	for i in 1:k
		tmp = vacuumstate(lattice)
		pos1 = index(lattice, i, band=band)
		for j in 1:k
			pos2 = index(lattice, j, band=band)
			coef = exp(corr[i, j]) - 1 
			if pos1==pos2
				t = ExpNTerm(pos1, coeff=coef)
			else
				t = ExpNTerm(pos1, pos2, coeff=coef)
			end
			apply!(t, tmp)
			canonicalize!(tmp, alg=Orthogonalize(TK.SVD(), trunc))
		end
		# println("-bond dimension of partial IF is ", bond_dimension(tmp))
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end

function hybriddynamics_2band!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	@assert lattice.bands == 2
	_hybriddynamics_1band!(gmps, lattice, corr1, 2, trunc=trunc)

	corr = corr1.data
	k = lattice.N
	for i in 1:k
		pos1 = index(lattice, i, band=1)
		tmp = vacuumstate(lattice)
		for j in 1:k, b2 in 1:lattice.bands
			coef = ifelse(b2==1, exp(corr[i, j]) - 1, exp(corr[i, j] + corr[j, i]) - 1)
			pos2 = index(lattice, j, band=b2)
			if pos1==pos2
				t = ExpNTerm(pos1, coeff=coef)
			else
				t = ExpNTerm(pos1, pos2, coeff=coef)
			end
			apply!(t, tmp)
			canonicalize!(tmp, alg=Orthogonalize(TK.SVD(), trunc))			
		end
		# println("--bond dimension of partial IF is ", bond_dimension(tmp))
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps	
end

# function _hybriddynamics_1band!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr1::ImagCorrelationFunction, band::Int; trunc::TruncationScheme=DefaultITruncation)
# 	corr = corr1.data
# 	k = lattice.N
# 	for i in 1:k
# 		coefs = [exp(corr[i, j]) - 1 for j in 1:k]
# 		row = index(lattice, i, band=band)
# 		cols = [index(lattice, j, band=band) for j in 1:k]
# 		tmp = partialdensempo(row, cols, coefs) * vacuumstate(lattice)
# 		mult!(gmps, tmp, trunc=trunc)
# 	end
# 	return gmps
# end

# function hybriddynamics_2band!(gmps::FockMPS, lattice::ImagFockLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
# 	@assert lattice.bands == 2
# 	_hybriddynamics_1band!(gmps, lattice, corr1, 2, trunc=trunc)

# 	corr = corr1.data
# 	k = lattice.N
# 	for i in 1:k
# 		pos1 = index(lattice, i, band=1)
# 		cols = Int[]
# 		coefs = scalartype(lattice)[]
# 		for j in 1:k, b2 in 1:lattice.bands
# 			coef = ifelse(b2==1, exp(corr[i, j]) - 1, exp(corr[i, j] + corr[j, i]) - 1)
# 			pos2 = index(lattice, j, band=b2)
# 			push!(cols, pos2)
# 			push!(coefs, coef)
# 		end
# 		tmp = partialdensempo(row, cols, coefs) * vacuumstate(lattice)
# 		mult!(gmps, tmp, trunc=trunc)
# 	end
# 	return gmps	
# end

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
			t = ExpNTerm(pos1, coeff=coef)
			apply!(t, gmps)
		else
			t = ExpNTerm(pos1, pos2, coeff=coef)
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
		t = ExpNTerm(pos1, pos2, coeff=coef)
		apply!(t, gmps)
		canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
	end
	return gmps	
end


