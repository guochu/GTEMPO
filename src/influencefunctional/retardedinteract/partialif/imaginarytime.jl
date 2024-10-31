"""
	retardedinteractdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, corr::ImagCorrelationFunction; band, trunc)

imaginary-time MPS-IF for a single band 
"""
function retardedinteractdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	corr = corr1.data
	k = lattice.k-1
	for i in 1:k
		tmp = partialif_retardedinteract(lattice, i, view(corr, i, 1:k), band=band)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps	
end

function partialif_retardedinteract(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1, trunc::TruncationScheme=DefaultMPOTruncation)
	if LayoutStyle(lattice) isa A1B1B1A1
		return _partialif_retardedinteract(lattice, i, cols; band=band, trunc=trunc)
	else
		lattice2 = similar(lattice, ordering=A1B1B1A1())
		gmps = _partialif_retardedinteract(lattice2, i, cols; band=band, trunc=trunc)
		perm = matchindices2(lattice, lattice2)
		return permute(gmps, perm, trunc=trunc)
	end
end

function _partialif_retardedinteract(lattice::ImagGrassmannLattice1Order{A2A2A1A1B2B2B1B1}, i::Int, cols::AbstractVector; band::Int=1, trunc::TruncationScheme=DefaultMPOTruncation)
	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	row = (index(lattice, i+1, band=band, conj=true), index(lattice, i, band=band, conj=false))
	# col_pos = [(index(lattice, j+1, band=band, conj=true), index(lattice, j, band=band, conj=false)) for j in length(cols):-1:1]
	col_pos = Tuple{Int, Int}[]
	cols2 = eltype(cols)[]
	for (j, c) in zip(length(cols):-1:1, reverse(cols))
		for b in 1:lattice.bands
			pos1, pos2 = index(lattice, j+1, band=b, conj=true), index(lattice, j, band=b, conj=false)
			push!(col_pos, (pos1, pos2))
			push!(cols2, c)
		end
	end
	p = sortperm(col_pos)
	mpo = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc)
	return mpo * vacuumstate(lattice)
end

function _partialif_retardedinteract(lattice::ImagGrassmannLattice1Order{A1B1B1A1}, i::Int, cols::AbstractVector; band::Int=1, trunc::TruncationScheme=DefaultMPOTruncation)
	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	row = (index(lattice, i+1, band=band, conj=true), index(lattice, i, band=band, conj=false))
	# col_pos = [(index(lattice, j+1, band=band, conj=true), index(lattice, j, band=band, conj=false)) for j in length(cols):-1:1]
	rev_cols = eltype(cols)[]
	for j in union(length(cols):-1:i+1, i-1:-1:1)
		push!(rev_cols, cols[j])
	end

	mps = vacuumstate(lattice)
	for band2 in 1:lattice.bands
		col_pos = Tuple{Int, Int}[]
		for j in union(length(cols):-1:i+1, i-1:-1:1)
			push!(col_pos, (index(lattice, j+1, band=band2, conj=true), index(lattice, j, band=band2, conj=false)))
		end		
		tmp = partialmpo_retardedinteract(row, col_pos, rev_cols, trunc=trunc) * vacuumstate(lattice)

		mult!(mps, tmp, trunc=trunc)
	end
	for band2 in 1:lattice.bands
		if band2 == band
			pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
			t = exp(GTerm(pos1, pos2, coeff=cols[i]))
		else
			pos2a, pos2b = index(lattice, i+1, conj=true, band=band2), index(lattice, i, conj=false, band=band2)
			t = exp(GTerm(row[1], row[2], pos2a, pos2b, coeff=cols[i]))
		end
		apply!(t, mps)
		canonicalize!(mps, alg=Orthogonalize(TK.SVD(), trunc))
	end
	return mps
end
