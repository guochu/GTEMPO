function partialif_retardedinteract(lattice::RealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; 
									b1::Symbol, band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	if LayoutStyle(lattice) isa A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1
		return _partialif_retardedinteract(lattice, i, cols_f, cols_b; b1=b1, band=band, trunc=trunc)
	else
		lattice2 = similar(lattice, ordering=A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1())
		gmps = _partialif_retardedinteract(lattice2, i, cols_f, cols_b; b1=b1, band=band, trunc=trunc)
		perm = matchindices2(lattice, lattice2)
		return permute(gmps, perm, trunc=trunc)
	end
end

function _partialif_retardedinteract(lattice::RealGrassmannLattice{A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1}, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; 
									b1::Symbol, band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	row = (index(lattice, i+1, band=band, conj=true, branch=b1), index(lattice, i, band=band, conj=false, branch=b1))
	# col_pos = [(index(lattice, j+1, band=band, conj=true), index(lattice, j, band=band, conj=false)) for j in length(cols):-1:1]
	col_pos = Tuple{Int, Int}[]
	cols2 = eltype(cols_f)[]
	for j in length(cols_f):-1:1
		for band2 in 1:lattice.bands
			for b2 in (:+, :-)
				c = (b2 == :+) ? cols_f[j] : cols_b[j]

				pos1, pos2 = index(lattice, j+1, band=band2, conj=true, branch=b2), index(lattice, j, band=band2, conj=false, branch=b2)
				push!(col_pos, (pos1, pos2))
				push!(cols2, c)

			end
		end
	end
	p = sortperm(col_pos)
	mpo = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc)
	return mpo * vacuumstate(lattice)
end

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