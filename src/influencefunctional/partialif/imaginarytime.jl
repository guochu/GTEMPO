"""
	hybriddynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, corr::ImagCorrelationFunction; band, trunc)

imaginary-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	corr = corr1.data
	k = lattice.k-1
	for i in 1:k
		tmp = partialinfluencefunctional(lattice, i+1, [0; view(corr, i, 1:k)], band=band)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end



function partialinfluencefunctional(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1)
	row = index(lattice, i, band=band, conj=true)
	col_pos = [index(lattice, j, band=band, conj=false) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end


"""
	retardedinteractdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, corr::ImagCorrelationFunction; band, trunc)

imaginary-time MPS-IF for a single band 
"""
function retardedinteractdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	
end