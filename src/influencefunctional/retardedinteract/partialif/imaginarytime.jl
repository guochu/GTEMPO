"""
	retardedinteractdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, corr::ImagCorrelationFunction; trunc)

imaginary-time MPS-IF for a single band 
"""
function retardedinteractdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	corr = corr1.data
	if lattice.bands == 1
		return retardedinteractdynamics_1band!(gmps, lattice, corr1, trunc=trunc)
	else
		if LayoutStyle(lattice) isa BandLocalLayout
			return retardedinteractdynamics_2band!(gmps, lattice, corr1, trunc=trunc)
		else

			k = lattice.k-1
			# 11 and 22
			for band in 1:lattice.bands
				for i in 1:k
					row = (index(lattice, i+1, band=band, conj=true), index(lattice, i, band=band, conj=false))
					col_pos = Tuple{Int, Int}[]
					cols = view(corr, i, 1:k)
					cols2 = eltype(cols)[]
					for (j, c) in zip(k:-1:1, reverse(cols))
						pos1, pos2 = index(lattice, j+1, band=band, conj=true), index(lattice, j, band=band, conj=false)
						push!(col_pos, (pos1, pos2))
						push!(cols2, c)
					end
					p = sortperm(col_pos)
					tmp = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
					mult!(gmps, tmp, trunc=trunc)
				end
			end
			# 12, naive approach
			for i in 1:k, j in 1:k
				pos1a, pos1b = index(lattice, i+1, conj=true, band=1), index(lattice, i, conj=false, band=1)
				pos2a, pos2b = index(lattice, j+1, conj=true, band=2), index(lattice, j, conj=false, band=2)
				coef = corr[i, j] + corr[j, i]
				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
				apply!(t, gmps)
				canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
			end
		end
	end
	return gmps	
end


function retardedinteractdynamics_1band!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	@assert lattice.bands == 1
	corr = corr1.data
	k = lattice.k-1
	band = 1
	for i in 1:k
		row = (index(lattice, i+1, band=band, conj=true), index(lattice, i, band=band, conj=false))
		col_pos = Tuple{Int, Int}[]
		cols = view(corr, i, 1:k)
		cols2 = eltype(cols)[]
		for (j, c) in zip(k:-1:1, reverse(cols))
			pos1, pos2 = index(lattice, j+1, band=band, conj=true), index(lattice, j, band=band, conj=false)
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
function retardedinteractdynamics_2band!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultMPOTruncation)
	@assert lattice.bands == 2
	(LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	corr = corr1.data
	k = lattice.k-1
	# IF 11 and IF 12
	band = 1
	for i in 1:k
		row = (index(lattice, i+1, band=band, conj=true), index(lattice, i, band=band, conj=false))
		col_pos = Tuple{Int, Int}[]
		cols2 = eltype(corr)[]
		for j in k:-1:1
			for b in 1:lattice.bands
				pos1, pos2 = index(lattice, j+1, band=b, conj=true), index(lattice, j, band=b, conj=false)
				push!(col_pos, (pos1, pos2))
				c = ifelse(b==band, corr[i, j], corr[i, j] + corr[j, i]) 
				push!(cols2, c)
			end
		end
		p = sortperm(col_pos)
		tmp = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
		mult!(gmps, tmp, trunc=trunc)
	end
	band = 2
	for i in 1:k
		row = (index(lattice, i+1, band=band, conj=true), index(lattice, i, band=band, conj=false))
		col_pos = Tuple{Int, Int}[]
		cols2 = eltype(corr)[]
		for j in k:-1:1
			pos1, pos2 = index(lattice, j+1, band=band, conj=true), index(lattice, j, band=band, conj=false)
			push!(col_pos, (pos1, pos2))
			push!(cols2, corr[i, j])
		end
		p = sortperm(col_pos)
		tmp = partialmpo_retardedinteract(row, col_pos[p], cols2[p], trunc=trunc) * vacuumstate(lattice)
		mult!(gmps, tmp, trunc=trunc)
	end

	return gmps
end