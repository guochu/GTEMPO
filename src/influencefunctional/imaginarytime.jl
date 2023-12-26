# partialinfluencefunctional(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; kwargs...) = pifrow(lattice, i, cols; kwargs...)
function partialinfluencefunctional2(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1)
	row = index(lattice, i, band=band, conj=true)
	col_pos = [index(lattice, j, band=band, conj=false) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end

# function pifrow(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1, trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	@assert (length(cols) <= lattice.k) && (i <= lattice.k)

# 	t = grassmanncreation()
# 	k = length(cols)
# 	data = Vector{typeof(t)}(undef, k)
# 	positions = Vector{Int}(undef, k)
# 	for (pos, j) in enumerate(k:-1:1)
# 		coef = (j < i) ? cols[j] : -cols[j]
# 		data[pos] = exp(t * coef)
# 		positions[pos] = index(lattice, j, band=band, conj=false)
# 	end
# 	positions = insert!(positions, k-i+2, index(lattice, i, band=band, conj=true))
# 	# println("1**********, positions ", positions)

# 	data = pif_util(data, k-i+1, trunc=DefaultMPOTruncation)
# 	mpo = QTerm(positions, data)
# 	mps = GrassmannMPS(length(lattice))
# 	out = mpo * mps
# 	return canonicalize!(out, alg=Orthogonalize(SVD(), trunc))
# end