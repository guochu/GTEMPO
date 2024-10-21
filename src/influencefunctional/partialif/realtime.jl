"""
	hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction; band::Int, trunc)

real-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	k = lattice.k
	for i in 1:k
		tmp1 = partialif_hybrid(lattice, i, view(η⁺⁺, i, 1:k), view(η⁺⁻, i, 1:k), b1=:+, band=band)
		tmp3 = partialif_hybrid(lattice, i, view(η⁻⁺, i, 1:k), view(η⁻⁻, i, 1:k), b1=:-, band=band)

		gmps = mult!(gmps, tmp1, trunc=trunc)
		gmps = mult!(gmps, tmp3, trunc=trunc)
	end
	return gmps		
end

function partialif_hybrid(lattice::RealGrassmannLattice, i::Int, cols::AbstractVector; b1::Symbol, b2::Symbol, band::Int=1)
	row = index(lattice, i, band=band, conj=true, branch=b1)
	col_pos = [index(lattice, j, band=band, conj=false, branch=b2) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end
function partialif_hybrid(lattice::RealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; b1::Symbol, band::Int=1)
	row = index(lattice, i, band=band, conj=true, branch=b1)
	cols = eltype(cols_f)[]
	col_pos = Int[]
	for j in length(cols_f):-1:1
		pos1, pos2 = index(lattice, j, band=band, conj=false, branch=:+), index(lattice, j, band=band, conj=false, branch=:-)
		push!(col_pos, pos1)
		push!(col_pos, pos2)
		push!(cols, cols_f[j])
		push!(cols, cols_b[j])
	end
	mpo = partialmpo(row, col_pos, cols)
	return mpo * vacuumstate(lattice)
end
function partialif_hybrid(lattice::RealGrassmannLattice, rows::AbstractVector, j::Int; b1::Symbol, b2::Symbol, band::Int=1)
	col = index(lattice, j, band=band, conj=false, branch=b2)
	row_pos = [index(lattice, i, band=band, conj=true, branch=b1) for i in length(rows):-1:1]
	mpo = partialmpo(col, row_pos, -reverse(rows))
	return mpo * vacuumstate(lattice)
end
function partialif_hybrid(lattice::RealGrassmannLattice, rows_f::AbstractVector, rows_b::AbstractVector, j::Int; b2::Symbol, band::Int=1)
	col = index(lattice, j, band=band, conj=false, branch=b2)
	rows = eltype(rows_f)[]
	row_pos = Int[]
	for i in length(rows_f):-1:1
		pos1, pos2 = index(lattice, i, band=band, conj=true, branch=:+), index(lattice, i, band=band, conj=true, branch=:-)
		push!(row_pos, pos1)
		push!(row_pos, pos2)
		push!(rows, -rows_f[i])
		push!(rows, -rows_b[i])
	end
	mpo = partialmpo(col, row_pos, rows)
	return mpo * vacuumstate(lattice)
end

function partialif_retardedinteract(lattice::RealGrassmannLattice{A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1}, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; 
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