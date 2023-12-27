function partialinfluencefunctional(lattice::RealGrassmannLattice, i::Int, cols::AbstractVector; fi::Bool, fj::Bool, band::Int=1)
	row = index(lattice, i, band=band, conj=true, forward=fi)
	col_pos = [index(lattice, j, band=band, conj=false, forward=fj) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end
function partialinfluencefunctional(lattice::RealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; fi::Bool, band::Int=1)
	row = index(lattice, i, band=band, conj=true, forward=fi)
	cols = eltype(cols_f)[]
	col_pos = Int[]
	for j in length(cols_f):-1:1
		pos1, pos2 = index(lattice, j, band=band, conj=false, forward=true), index(lattice, j, band=band, conj=false, forward=false)
		push!(col_pos, pos1)
		push!(col_pos, pos2)
		push!(cols, cols_f[j])
		push!(cols, cols_b[j])
	end
	mpo = partialmpo(row, col_pos, cols)
	return mpo * vacuumstate(lattice)
end
function partialinfluencefunctional(lattice::RealGrassmannLattice, rows::AbstractVector, j::Int; fi::Bool, fj::Bool, band::Int=1)
	col = index(lattice, j, band=band, conj=false, forward=fj)
	row_pos = [index(lattice, i, band=band, conj=true, forward=fi) for i in length(rows):-1:1]
	mpo = partialmpo(col, row_pos, -reverse(rows))
	return mpo * vacuumstate(lattice)
end
function partialinfluencefunctional(lattice::RealGrassmannLattice, rows_f::AbstractVector, rows_b::AbstractVector, j::Int; fj::Bool, band::Int=1)
	col = index(lattice, j, band=band, conj=false, forward=fj)
	rows = eltype(rows_f)[]
	row_pos = Int[]
	for i in length(rows_f):-1:1
		pos1, pos2 = index(lattice, i, band=band, conj=true, forward=true), index(lattice, i, band=band, conj=true, forward=false)
		push!(row_pos, pos1)
		push!(row_pos, pos2)
		push!(rows, -rows_f[i])
		push!(rows, -rows_b[i])
	end
	mpo = partialmpo(col, row_pos, rows)
	return mpo * vacuumstate(lattice)
end

