


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
			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
		end
	end
	return gmps
end

function retardedinteractdynamics_2band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; 
											trunc::TruncationScheme=DefaultMPOTruncation)
	@assert lattice.bands == 2
	for band in 1:lattice.bands
		_retardedinteractdynamics_1band!(gmps, lattice, corr, band, trunc=trunc)
	end

	for i in 1:lattice.kt-1, b1 in (:+, :-)
		pos1a, pos1b, c1 = get_pair_pos(lattice, i, 1, b1)
		for j in 1:lattice.kt-1, b2 in (:+, :-)
			pos2a, pos2b, c2 = get_pair_pos(lattice, j, 2, b2)
			c = 2 * index(corr, i, j, b1=b1, b2=b2)  * c1 * c2
			t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=c))
			apply!(t, gmps)
			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
		end
	end
	return gmps
end


# function retardedinteractdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
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


# function _retardedinteractdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
# 	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
# 	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
# 	if lattice.bands == 1
# 		return retardedinteractdynamics_1band!(gmps, lattice, corr, trunc=trunc)
# 	else

# 		return retardedinteractdynamics_2band!(gmps, lattice, corr, trunc=trunc)
# 	end
# end

# function retardedinteractdynamics_1band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; 
# 											trunc::TruncationScheme=DefaultITruncation)	
# 	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
# 	@assert lattice.bands == 1
# 	k = lattice.k-1
# 	band = 1
# 	for i in 1:k, b1 in (:+, :-)
# 		row1, row2, c1 = get_pair_pos(lattice, i, band, b1) 

# 		col_pos = Tuple{Int, Int}[]
# 		cols2 = scalartype(lattice)[]
# 		for j in k:-1:1, b2 in (:+, :-)
# 			pos1, pos2, c2 = get_pair_pos(lattice, j, band, b2)
# 			c = index(corr, i, j, b1=b1, b2=b2) * c1 *c2
# 			push!(col_pos, (pos1, pos2))
# 			push!(cols2, c)
# 		end
# 		p = sortperm(col_pos)
# 		println(col_pos[p])
# 		tmp = partialmpo_retardedinteract((row1, row2), col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
# 		mult!(gmps, tmp, trunc=trunc)
# 	end
# 	return gmps
# end

# # only applicable for BandLocalLayout
# function retardedinteractdynamics_2band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultMPOTruncation)
# 	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
# 	@assert lattice.bands == 2
# 	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
# 	k = lattice.k-1
# 	# IF 11 and IF 12
# 	band = 1
# 	for i in 1:k, b1 in (:+, :-)
# 		row1, row2, c1 = get_pair_pos(lattice, i, band, b1) 
# 		col_pos = Tuple{Int, Int}[]
# 		cols2 = scalartype(lattice)[]
# 		for j in k:-1:1, b in 1:lattice.bands, b2 in (:+, :-)
# 			pos1, pos2, c2 = get_pair_pos(lattice, j, b, b2)
# 			push!(col_pos, (pos1, pos2))
# 			c = ifelse(b==band, index(corr, i, j, b1=b1, b2=b2), index(corr, i, j, b1=b1, b2=b2) + index(corr, j, i, b1=b1, b2=b2)) 
# 			push!(cols2, c*c1*c2)
# 		end
# 		p = sortperm(col_pos)
# 		tmp = partialmpo_retardedinteract((row1, row2), col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
# 		mult!(gmps, tmp, trunc=trunc)
# 	end
# 	band = 2
# 	for i in 1:k, b1 in (:+, :-)
# 		row1, row2, c1 = get_pair_pos(lattice, i, band, b1)
# 		col_pos = Tuple{Int, Int}[]
# 		cols2 = scalartype(lattice)[]
# 		for j in k:-1:1, b2 in (:+, :-)
# 			pos1, pos2, c2 = get_pair_pos(lattice, j, band, b2)
# 			push!(col_pos, (pos1, pos2))
# 			push!(cols2, index(corr, i, j, b1=b1, b2=b2)*c1*c2)
# 		end
# 		p = sortperm(col_pos)
# 		tmp = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
# 		mult!(gmps, tmp, trunc=trunc)
# 	end

# 	return gmps
# end



function get_pair_pos(lattice::AbstractGrassmannLattice, j::Int, band::Int, b::Symbol)
	if b == :-
		return index(lattice, j+1, band=band, conj=false, branch=b), index(lattice, j, band=band, conj=true, branch=b), -1
	else
		return index(lattice, j+1, band=band, conj=true, branch=b), index(lattice, j, band=band, conj=false, branch=b), 1
	end
end

