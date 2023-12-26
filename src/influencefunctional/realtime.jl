# partialinfluencefunctional(lattice::RealGrassmannLattice, i::Int, cols::AbstractVector; kwargs...) = pifrow(lattice, i, cols; kwargs...)
# partialinfluencefunctional(lattice::RealGrassmannLattice, rows::AbstractVector, j::Int; kwargs...) = pifcol(lattice, rows, j; kwargs...)
function partialinfluencefunctional2(lattice::RealGrassmannLattice, i::Int, cols::AbstractVector; fi::Bool, fj::Bool, band::Int=1)
	row = index(lattice, i, band=band, conj=true, forward=fi)
	col_pos = [index(lattice, j, band=band, conj=false, forward=fj) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end
function partialinfluencefunctional2(lattice::RealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; fi::Bool, band::Int=1)
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
function partialinfluencefunctional2(lattice::RealGrassmannLattice, rows::AbstractVector, j::Int; fi::Bool, fj::Bool, band::Int=1)
	col = index(lattice, j, band=band, conj=false, forward=fj)
	row_pos = [index(lattice, i, band=band, conj=true, forward=fi) for i in length(rows):-1:1]
	mpo = partialmpo(col, row_pos, -reverse(rows))
	return mpo * vacuumstate(lattice)
end
# function partialinfluencefunctional2(lattice::RealGrassmannLattice, rows_f::AbstractVector, rows_b::AbstractVector, j::Int; fj::Bool, band::Int=1)
# 	@assert (lattice.ordering in (ABBA(), AABB()))
# 	col = index(lattice, j, band=band, conj=false, forward=fj)
# 	rows = eltype(rows_f)[]
# 	row_pos = Int[]
# 	for i in length(rows_f):-1:1
# 		if lattice.ordering isa AABB
# 			pos1, pos2 = index(lattice, i, band=band, conj=true, forward=true), index(lattice, i, band=band, conj=true, forward=false)
# 			push!(row_pos, pos1)
# 			push!(row_pos, pos2)
# 			push!(rows, -rows_f[i])
# 			push!(rows, -rows_b[i])
# 		elseif lattice.ordering isa ABBA
# 			pos1, pos2 = index(lattice, i, band=band, conj=true, forward=true), index(lattice, i, band=band, conj=true, forward=false)
# 			push!(row_pos, pos2)
# 			push!(row_pos, pos1)
# 			push!(rows, -rows_b[i])
# 			push!(rows, -rows_f[i])
# 		else			
# 			error("partialinfluencefunctional2 not implemented for GrassmannOrdering type $(typeof(lattice.ordering))")
# 		end
# 	end
# 	# println(row_pos)
# 	mpo = partialmpo(col, row_pos, rows)
# 	return mpo * vacuumstate(lattice)
# end
function partialinfluencefunctional2(lattice::RealGrassmannLattice, rows_f::AbstractVector, rows_b::AbstractVector, j::Int; fj::Bool, band::Int=1)
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
	# println(row_pos)
	mpo = partialmpo(col, row_pos, rows)
	return mpo * vacuumstate(lattice)
end

# """
# 	pifcol(lattice::RealGrassmannLattice2Order, i::Int, cols::AbstractVector; kwargs...)

# The influence functional IF[rows, j] as MPO
# """
# function pifcol(lattice::RealGrassmannLattice, rows::AbstractVector, j::Int; 
# 										fi::Bool, fj::Bool, band::Int=1, trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	@assert (length(rows) <= lattice.k) && (j <= lattice.k)
# 	t = complex(grassmanncreation())
# 	k = length(rows)

# 	# all the rows, the k-th column
# 	data = Vector{typeof(t)}(undef, k)
# 	positions = Vector{Int}(undef, k)
# 	pos_j = index(lattice, j, band=band, conj=false, forward=fj)
# 	for (pos, i) in enumerate(k:-1:1)
# 		pos_i = index(lattice, i, band=band, conj=true, forward=fi)
# 		# coef = (j < i) ? rows[i] : -rows[i]
# 		coef = (pos_j > pos_i) ? rows[i] : -rows[i]
# 		data[pos] = exp(t * coef)
# 		positions[pos] = pos_i
# 	end
# 	insert_pos = findfirst(x->x>pos_j, positions)
# 	if isnothing(insert_pos)
# 		insert_pos = k+1
# 	end
# 	positions = insert!(positions, insert_pos, pos_j)
# 	data = pif_util(data, insert_pos-1, trunc=DefaultMPOTruncation)
# 	mpo = QTerm(positions, data)
# 	mps = GrassmannMPS(ComplexF64, length(lattice))
# 	out = mpo * mps	

# 	return canonicalize!(out, alg=Orthogonalize(SVD(), trunc))
# end

# """
# 	pifrow(lattice::RealGrassmannLattice2Order, i::Int, cols::AbstractVector; kwargs...)

# The influence functional IF[i, cols] as MPO
# """
# function pifrow(lattice::RealGrassmannLattice, i::Int, cols::AbstractVector; 
# 										fi::Bool, fj::Bool, band::Int=1, trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	@assert (length(cols) <= lattice.k) && (i <= lattice.k)
# 	t = complex(grassmanncreation())
# 	k = length(cols)

# 	# the k-th row, all columns
# 	data = Vector{typeof(t)}(undef, k)
# 	positions = Vector{Int}(undef, k)
# 	pos_i = index(lattice, i, band=band, conj=true, forward=fi)
# 	for (pos, j) in enumerate(k:-1:1)
# 		pos_j = index(lattice, j, band=band, conj=false, forward=fj)
# 		# coef = (j < i) ? cols[j] : -cols[j]
# 		coef = (pos_j > pos_i) ? cols[j] : -cols[j]
# 		data[pos] = exp(t * coef)
# 		positions[pos] = pos_j
# 	end
# 	insert_pos = findfirst(x->x>pos_i, positions)
# 	if isnothing(insert_pos)
# 		insert_pos = k+1
# 	end
# 	positions = insert!(positions, insert_pos, pos_i)
# 	data = pif_util(data, insert_pos-1, trunc=DefaultMPOTruncation)
# 	mpo = QTerm(positions, data)
# 	mps = GrassmannMPS(ComplexF64, length(lattice))
# 	out = mpo * mps	
# 	return canonicalize!(out, alg=Orthogonalize(SVD(), trunc))
# end

