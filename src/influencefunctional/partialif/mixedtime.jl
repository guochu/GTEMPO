"""
	hybriddynamics(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::MixedCorrelationFunction; band::Int, trunc)

real-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::AbstractMixedCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	kt, Nτ = lattice.kt, lattice.Nτ
	for b1 in (:+, :-)
		for i in 1:kt
			cols_f = [index(corr, i, j, b1=b1, b2=:+) for j in 1:kt]
			cols_b = [index(corr, i, j, b1=b1, b2=:-) for j in 1:kt]
			cols_i = [index(corr, i, j, b1=b1, b2=:τ) for j in 1:Nτ]
			tmp = partialinfluencefunctional(lattice, i, cols_f, cols_b, cols_i, b1=b1, band=band)
			gmps = mult!(gmps, tmp, trunc=trunc)
		end
	end
	b1 = :τ
	for i in 1:Nτ
		cols_f = [index(corr, i, j, b1=b1, b2=:+) for j in 1:kt]
		cols_b = [index(corr, i, j, b1=b1, b2=:-) for j in 1:kt]
		cols_i = [index(corr, i, j, b1=b1, b2=:τ) for j in 1:Nτ]
		tmp = partialinfluencefunctional(lattice, i, cols_f, cols_b, cols_i, b1=b1, band=band)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps		
end

function partialinfluencefunctional(lattice::MixedGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector, cols_i::AbstractVector; b1::Symbol, band::Int=1)
	@assert length(cols_f) == length(cols_b)
	@assert length(cols_i) == lattice.Nτ
	if b1 == :τ
		i = i + 1
	end
	row = index(lattice, i, band=band, conj=true, branch=b1)
	cols = eltype(cols_f)[]
	col_pos = Int[]

	for j in 1:length(cols_i)
		pos = index(lattice, j+1, band=band, conj=false, branch=:τ)
		push!(col_pos, pos)
		push!(cols, cols_i[j])
	end

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

# function partialinfluencefunctional(lattice::MixedGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector, cols_i::AbstractVector; b1::Symbol, band::Int=1)
# 	@assert length(cols_f) == length(cols_b)
# 	@assert length(cols_i) == lattice.Nτ
# 	# if b1 == :τ
# 	# 	i = i + 1
# 	# end
# 	row = index(lattice, i, band=band, conj=true, branch=b1)
# 	cols = eltype(cols_f)[]
# 	col_pos = Int[]

# 	for j in 1:length(cols_i)
# 		pos = index(lattice, j, band=band, conj=false, branch=:τ)
# 		push!(col_pos, pos)
# 		push!(cols, cols_i[j])
# 	end

# 	for j in length(cols_f):-1:1
# 		pos1, pos2 = index(lattice, j, band=band, conj=false, branch=:+), index(lattice, j, band=band, conj=false, branch=:-)
# 		push!(col_pos, pos1)
# 		push!(col_pos, pos2)
# 		push!(cols, cols_f[j])
# 		push!(cols, cols_b[j])
# 	end
	
# 	mpo = partialmpo(row, col_pos, cols)
# 	return mpo * vacuumstate(lattice)
# end