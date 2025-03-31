"""
	focktograssmann(lattice::AbstractFockLattice, mps::FockMPS; trunc::TruncationScheme=DefaultKTruncation)
"""
function focktograssmann(lattice::AbstractFockLattice, mps::FockMPS)
	lattice2 = similargrassmannlattice(lattice)
	mps2 = vacuumstate(lattice2)
	for band in 1:lattice.bands, b in branches(lattice)
		if b == :τ
			k = lattice.Nτ
		else
			k = lattice.Nt
		end
		for j in 1:k
			focktograssmann_one_step!(lattice2, mps2, lattice, mps, j, band, b)
		end
	end
	new_scale = scaling(mps)^(length(mps) / length(mps2))
	setscaling!(mps2, new_scale)
	return lattice2, mps2
end


function focktograssmann(::Type{O}, lattice::AbstractFockLattice, mps::FockMPS) where {O<:GrassmannOrdering}
	lattice2, mps2 = focktograssmann(lattice, mps)
	if !(lattice2.ordering isa O)
		lattice2, mps2 = changeordering(O, lattice2, mps2)
	end
	return lattice2, mps2
end


function focktograssmann_one_step!(lattice2::AbstractGrassmannLattice, mps2::GrassmannMPS, lattice::AbstractFockLattice, mps::FockMPS, j::Int, band::Int, b::Symbol)
	ph = grassmannpspace()
	pos = index(lattice, j, band=band, branch=b)
	t = mps[pos]
	if b == :-
		pos1, pos2 = index(lattice2, j, band=band, branch=b, conj=true), index(lattice2, j+1, band=band, branch=b, conj=false)
	else
		pos1, pos2 = index(lattice2, j+1, band=band, branch=b, conj=true), index(lattice2, j, band=band, branch=b, conj=false)
	end
	# println("j=", j, ",band=", band, ",branch=", b, ",pos1=", pos1, ",pos2=", pos2)
	(pos1 + 1 == pos2) || error("something wrong")
	left = Rep[ℤ₂](0=>space_l(t))
	right = Rep[ℤ₂](0=>space_r(t))
	t2 = zeros(scalartype(mps), left ⊗ ph ← ph' ⊗ right)
	copy!(block(t2, Z2Irrep(0)), t[:, 1, :])
	copy!(block(t2, Z2Irrep(1)), t[:, 2, :])
	x, y = rightorth!(t2, alg=LQ())
	mps2[pos1] = x
	mps2[pos2] = permute(y, (1,2), (3,))
end