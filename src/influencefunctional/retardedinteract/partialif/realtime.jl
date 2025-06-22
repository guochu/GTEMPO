function retardedinteractdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	if lattice.bands == 1
		return retardedinteractdynamics_1band!(gmps, lattice, corr, trunc=trunc)
	else
		return retardedinteractdynamics_2band!(gmps, lattice, corr, trunc=trunc)
	end
end

function retardedinteractdynamics_1band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; kwargs...)
	@assert lattice.bands == 1
	return _retardedinteractdynamics_1band!(gmps, lattice, corr, 1; kwargs...)
end 


function _retardedinteractdynamics_1band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction, band::Int; 
											trunc::TruncationScheme=DefaultITruncation)	
	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	alg = Orthogonalize(TK.SVD(), trunc)
	for i in 1:lattice.kt-1, b1 in (:+, :-)
		tmp = vacuumstate(lattice)
		pos1a, pos1b, c1 = get_pair_pos(lattice, i, band, b1)
		for j in 1:lattice.kt-1, b2 in (:+, :-)
			pos2a, pos2b, c2 = get_pair_pos(lattice, j, band, b2)
			c = index(corr, i, j, b1=b1, b2=b2) * c1 * c2
			if (i == j) && (b1 == b2)
				t = exp(GTerm(pos1a, pos1b, coeff=c))
			else
				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=c))
			end
			apply!(t, tmp)
			canonicalize!(tmp, alg=alg)
		end
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end

function retardedinteractdynamics_2band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; 
											trunc::TruncationScheme=DefaultMPOTruncation)
	@assert lattice.bands == 2
	_retardedinteractdynamics_1band!(gmps, lattice, corr, 2, trunc=trunc)

	alg = Orthogonalize(TK.SVD(), trunc)
	for i in 1:lattice.kt-1, b1 in branches(lattice)
		tmp = vacuumstate(lattice)
		pos1a, pos1b, c1 = get_pair_pos(lattice, i, 1, b1)
		for j in 1:lattice.kt-1, b2 in branches(lattice), band2 in 1:lattice.bands
			pos2a, pos2b, c2 = get_pair_pos(lattice, j, band2, b2)
			if band2 == 1
				c = index(corr, i, j, b1=b1, b2=b2) * c1 * c2
			else
				c = (index(corr, i, j, b1=b1, b2=b2) + index(corr, j, i, b1=b2, b2=b1) )  * c1 * c2
			end
			if pos1a == pos2a
				@assert pos1b == pos2b
				t = exp(GTerm(pos1a, pos1b, coeff=c))
			else
				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=c))
			end
			apply!(t, tmp)
			canonicalize!(tmp, alg=alg)
		end
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end

# function retardedinteractdynamics_2band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; 
# 											trunc::TruncationScheme=DefaultMPOTruncation)
# 	@assert lattice.bands == 2
# 	_retardedinteractdynamics_1band!(gmps, lattice, corr, 2, trunc=trunc)

# 	alg = Orthogonalize(TK.SVD(), trunc)
# 	for i in 1:lattice.kt-1, b1 in (:+, :-)
# 		tmp = vacuumstate(lattice)
# 		pos1a, pos1b, c1 = get_pair_pos(lattice, i, 1, b1)
# 		for j in 1:lattice.kt-1, b2 in (:+, :-)
# 			pos2a, pos2b, c2 = get_pair_pos(lattice, j, 2, b2)
# 			c = (index(corr, i, j, b1=b1, b2=b2) + index(corr, j, i, b1=b2, b2=b1) )  * c1 * c2
# 			t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=c))
# 			apply!(t, tmp)
# 			canonicalize!(tmp, alg=alg)
# 		end
# 		for j in 1:lattice.kt-1, b2 in (:+, :-)
# 			pos2a, pos2b, c2 = get_pair_pos(lattice, j, 1, b2)
# 			c = index(corr, i, j, b1=b1, b2=b2) * c1 * c2
# 			if (i == j) && (b1 == b2)
# 				t = exp(GTerm(pos1a, pos1b, coeff=c))
# 			else
# 				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=c))
# 			end
# 			apply!(t, tmp)
# 			canonicalize!(tmp, alg=alg)
# 		end		
# 		mult!(gmps, tmp, trunc=trunc)
# 	end
# 	return gmps
# end

# naive implementation with N^2 gate operations
function retardedinteractdynamics_naive!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	if lattice.bands == 1
		return retardedinteractdynamics_1band_naive!(gmps, lattice, corr, trunc=trunc)
	else
		return retardedinteractdynamics_2band_naive!(gmps, lattice, corr, trunc=trunc)
	end
end

function retardedinteractdynamics_1band_naive!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; kwargs...)
	@assert lattice.bands == 1
	return _retardedinteractdynamics_1band_naive!(gmps, lattice, corr, 1; kwargs...)
end 


function _retardedinteractdynamics_1band_naive!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction, band::Int; 
												trunc::TruncationScheme=DefaultITruncation)	
	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	alg = Orthogonalize(TK.SVD(), trunc)
	for i in 1:lattice.kt-1, b1 in (:+, :-)
		pos1a, pos1b, c1 = get_pair_pos(lattice, i, band, b1)
		for j in 1:lattice.kt-1, b2 in (:+, :-)
			pos2a, pos2b, c2 = get_pair_pos(lattice, j, band, b2)
			c = index(corr, i, j, b1=b1, b2=b2) * c1 * c2
			if (i == j) && (b1 == b2)
				t = exp(GTerm(pos1a, pos1b, coeff=c))
			else
				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=c))
			end
			apply!(t, gmps)
			canonicalize!(gmps, alg=alg)
		end
	end
	return gmps
end

function retardedinteractdynamics_2band_naive!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; 
											trunc::TruncationScheme=DefaultMPOTruncation)
	@assert lattice.bands == 2
	for band in 1:lattice.bands
		_retardedinteractdynamics_1band_naive!(gmps, lattice, corr, band, trunc=trunc)
	end
	alg = Orthogonalize(TK.SVD(), trunc)
	for i in 1:lattice.kt-1, b1 in (:+, :-)
		pos1a, pos1b, c1 = get_pair_pos(lattice, i, 1, b1)
		for j in 1:lattice.kt-1, b2 in (:+, :-)
			pos2a, pos2b, c2 = get_pair_pos(lattice, j, 2, b2)
			c = (index(corr, i, j, b1=b1, b2=b2) + index(corr, j, i, b1=b2, b2=b1))  * c1 * c2
			t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=c))
			apply!(t, gmps)
			canonicalize!(gmps, alg=alg)
		end
	end
	return gmps
end



function get_pair_pos(lattice::AbstractGrassmannLattice, j::Int, band::Int, b::Symbol)
	if b == :-
		return index(lattice, j+1, band=band, conj=false, branch=b), index(lattice, j, band=band, conj=true, branch=b), -1
	else
		return index(lattice, j+1, band=band, conj=true, branch=b), index(lattice, j, band=band, conj=false, branch=b), 1
	end
end

