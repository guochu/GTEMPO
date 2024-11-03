
function retardedinteractdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	if LayoutStyle(lattice) isa BandLocalLayout
		return _retardedinteractdynamics!(gmps, lattice, corr, trunc=trunc)
	else
		lattice2 = similar(lattice, ordering=A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1())
		gmps2 = _retardedinteractdynamics!(vacuumstate(lattice2), lattice2, corr; trunc=trunc)
		perm = matchindices2(lattice, lattice2)
		gmps2 = permute(gmps2, perm, trunc=trunc)
		return mult!(gmps, gmps2, trunc=trunc)
	end
end


function _retardedinteractdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	if lattice.bands == 1
		return retardedinteractdynamics_1band!(gmps, lattice, corr, trunc=trunc)
	else

		return retardedinteractdynamics_2band!(gmps, lattice, corr, trunc=trunc)
	end
end

function retardedinteractdynamics_1band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; 
											trunc::TruncationScheme=DefaultITruncation)	
	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	@assert lattice.bands == 1
	k = lattice.k-1
	band = 1
	for i in 1:k, b1 in (:+, :-)
		row = (index(lattice, i+1, band=band, conj=true, branch=b1), index(lattice, i, band=band, conj=false, branch=b1))

		col_pos = Tuple{Int, Int}[]
		cols2 = scalartype(lattice)[]
		for j in k:-1:1, b2 in (:+, :-)
			pos1, pos2 = index(lattice, j+1, band=band, conj=true, branch=b2), index(lattice, j, band=band, conj=false, branch=b2)
			c = index(corr, i, j, b1=b1, b2=b2)
			push!(col_pos, (pos1, pos2))
			push!(cols2, c)
		end
		p = sortperm(col_pos)
		tmp = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
		mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end

# only applicable for BandLocalLayout
function retardedinteractdynamics_2band!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; trunc::TruncationScheme=DefaultMPOTruncation)
	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	@assert lattice.bands == 2
	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	k = lattice.k-1
	# IF 11 and IF 12
	band = 1
	for i in 1:k, b1 in (:+, :-)
		row = (index(lattice, i+1, band=band, conj=true, branch=b1), index(lattice, i, band=band, conj=false, branch=b1))
		col_pos = Tuple{Int, Int}[]
		cols2 = scalartype(lattice)[]
		for j in k:-1:1, b in 1:lattice.bands, b2 in (:+, :-)
			pos1, pos2 = index(lattice, j+1, band=b, conj=true, branch=b2), index(lattice, j, band=b, conj=false, branch=b2)
			push!(col_pos, (pos1, pos2))
			c = ifelse(b==band, index(corr, i, j, b1=b1, b2=b2), index(corr, i, j, b1=b1, b2=b2) + index(corr, j, i, b1=b1, b2=b2)) 
			push!(cols2, c)
		end
		p = sortperm(col_pos)
		tmp = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
		mult!(gmps, tmp, trunc=trunc)
	end
	band = 2
	for i in 1:k, b1 in (:+, :-)
		row = (index(lattice, i+1, band=band, conj=true, branch=b1), index(lattice, i, band=band, conj=false, branch=b1))
		col_pos = Tuple{Int, Int}[]
		cols2 = scalartype(lattice)[]
		for j in k:-1:1, b2 in (:+, :-)
			pos1, pos2 = index(lattice, j+1, band=band, conj=true, branch=b2), index(lattice, j, band=band, conj=false, branch=b2)
			push!(col_pos, (pos1, pos2))
			push!(cols2, index(corr, i, j, b1=b1, b2=b2))
		end
		p = sortperm(col_pos)
		tmp = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
		mult!(gmps, tmp, trunc=trunc)
	end

	return gmps
end


# function partialif_retardedinteract(lattice::RealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; 
# 									b1::Symbol, band::Int=1, trunc::TruncationScheme=DefaultITruncation)
# 	if LayoutStyle(lattice) isa A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1
# 		return _partialif_retardedinteract(lattice, i, cols_f, cols_b; b1=b1, band=band, trunc=trunc)
# 	else
# 		lattice2 = similar(lattice, ordering=A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1())
# 		gmps = _partialif_retardedinteract(lattice2, i, cols_f, cols_b; b1=b1, band=band, trunc=trunc)
# 		perm = matchindices2(lattice, lattice2)
# 		return permute(gmps, perm, trunc=trunc)
# 	end
# end

# function _partialif_retardedinteract(lattice::RealGrassmannLattice{A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1}, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; 
# 									b1::Symbol, band::Int=1, trunc::TruncationScheme=DefaultITruncation)
# 	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
# 	row = (index(lattice, i+1, band=band, conj=true, branch=b1), index(lattice, i, band=band, conj=false, branch=b1))
# 	# col_pos = [(index(lattice, j+1, band=band, conj=true), index(lattice, j, band=band, conj=false)) for j in length(cols):-1:1]
# 	col_pos = Tuple{Int, Int}[]
# 	cols2 = eltype(cols_f)[]
# 	for j in length(cols_f):-1:1
# 		for band2 in 1:lattice.bands
# 			for b2 in (:+, :-)
# 				c = (b2 == :+) ? cols_f[j] : cols_b[j]

# 				pos1, pos2 = index(lattice, j+1, band=band2, conj=true, branch=b2), index(lattice, j, band=band2, conj=false, branch=b2)
# 				push!(col_pos, (pos1, pos2))
# 				push!(cols2, c)

# 			end
# 		end
# 	end
# 	p = sortperm(col_pos)
# 	mpo = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc)
# 	return mpo * vacuumstate(lattice)
# end

# """
# 	partialif_retardedinteract(lattice::RealGrassmannLattice{A1a1B1b1b1B1a1A1}, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; b1::Symbol, band::Int=1)
# """
# function partialif_retardedinteract(lattice::RealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; 
# 									b1::Symbol, band::Int=1, trunc::TruncationScheme=DefaultITruncation)
# 	row = (index(lattice, i+1, band=band, conj=true, branch=b1), index(lattice, i, band=band, conj=false, branch=b1))
# 	rev_cols = eltype(cols_f)[]
# 	for j in union(length(cols_f):-1:i+1, i-1:-1:1)
# 		push!(rev_cols, cols_f[j])
# 		push!(rev_cols, cols_b[j])
# 	end

# 	mps = vacuumstate(lattice)
# 	for band2 in 1:lattice.bands
# 		col_pos = Tuple{Int, Int}[]
# 		for j in union(length(cols_f):-1:i+1, i-1:-1:1)
# 			for b in (:+, :-)
# 				push!(col_pos, (index(lattice, j+1, band=band2, conj=true, branch=b), index(lattice, j, band=band2, conj=false, branch=b)))
# 			end
# 		end	
# 		p = sortperm(col_pos)

# 		tmp = partialmpo_retardedinteract(row, col_pos[p], rev_cols[p], trunc=trunc) * vacuumstate(lattice)

# 		mult!(mps, tmp, trunc=trunc)
# 	end

# 	for band2 in 1:lattice.bands, b in (:+, :-)
# 		coef = (b == :+) ? cols_f[i] : cols_b[i]
# 		if (band2 == band) && (b1 == b)
# 			pos1, pos2 = index(lattice, i+1, conj=true, band=band, branch=b), index(lattice, i, conj=false, band=band, branch=b)
# 			t = exp(GTerm(pos1, pos2, coeff=coef))
# 		else
# 			pos2a, pos2b = index(lattice, i+1, conj=true, band=band2, branch=b), index(lattice, i, conj=false, band=band2, branch=b)
# 			t = exp(GTerm(row[1], row[2], pos2a, pos2b, coeff=coef))
# 		end
# 		apply!(t, mps)
# 		canonicalize!(mps, alg=Orthogonalize(TK.SVD(), trunc))
# 	end

# 	return mps
# end