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

function hybriddynamics_naive!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	k = lattice.k
	for i in 1:k
		tmp1 = partialif_hybrid_naive(lattice, i, view(η⁺⁺, i, 1:k), view(η⁺⁻, i, 1:k), b1=:+, band=band, trunc=trunc)
		tmp3 = partialif_hybrid_naive(lattice, i, view(η⁻⁺, i, 1:k), view(η⁻⁻, i, 1:k), b1=:-, band=band, trunc=trunc)

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
	@assert length(cols_f) == length(cols_b)
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
function partialif_hybrid_naive(lattice::RealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; 
							b1::Symbol, band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	@assert length(cols_f) == length(cols_b)
	alg = Orthogonalize(TK.SVD(), trunc)
	gmps = vacuumstate(lattice)
	for j in length(cols_f):-1:1
		pos1, pos2 = index(lattice, i, conj=true, branch=b1, band=band), index(lattice, j, conj=false, branch=:+, band=band)
		t = exp(GTerm(pos1, pos2, coeff=cols_f[j]))
		apply!(t, gmps)
		canonicalize!(gmps, alg=alg)
		pos1, pos2 = index(lattice, i, conj=true, branch=b1, band=band), index(lattice, j, conj=false, branch=:-, band=band)
		t = exp(GTerm(pos1, pos2, coeff=cols_b[j]))
		apply!(t, gmps)
		canonicalize!(gmps, alg=alg)
	end
	return gmps
end

# reverse col and row, only used in stepper
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