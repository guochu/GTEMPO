function partialinfluencefunctional(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1)
	row = index(lattice, i, band=band, conj=true)
	col_pos = [index(lattice, j, band=band, conj=false) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end
