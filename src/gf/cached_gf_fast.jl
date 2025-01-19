function cached_gf_fast(lattice::AbstractGrassmannLattice, N::Int, A::Union{GrassmannMPS, Vector}, Bs::GrassmannMPS...; 
					cache::GTEMPO.TwosideExpectationCache=environments(lattice, A, Bs...), 
					branch::Symbol, c1::Bool=true, c2::Bool=false, band::Int=1, kwargs...)
	# calculated G11
	a, b = ContourIndex(1, conj=c1, branch=branch, band=band), ContourIndex(1, conj=c2, branch=branch, band=band)
	g0 = cached_gf(lattice, a, b, A, Bs...; cache=cache, band=band, kwargs...)
	GFt = Vector{scalartype(lattice)}(undef, N)
	GFt[1] = g0
	(N == 1) && return GFt
	
	# calculated G12
	pos1 = index(lattice, 1, conj=c1, band=band, branch=branch) # left boundary position
	pos2 = index(lattice, 2, conj=c2, band=band, branch=branch) # right boundary position
	m = convert(PartialMPO, GTerm(pos1, pos2, coeff=1))
	A2 = m * cache.A

	j = pos2pairindex(positions(m)[1])
	k = pos2pairindex(positions(m)[end])
	left = leftenv(cache, j)
	right = rightenv(cache, k)
	for tj in j:k-1
		left = left * GrassmannTransferMatrix(tj, A2, cache.Bs...)
	end
	pos_j = k
	left′ = left * GrassmannTransferMatrix(k, A2, cache.Bs...)
	GFt[2] = contract_center(left′, right) / Zvalue(cache)

	# calculated G1k
	for i in 3:N+1
		pos2 = index(lattice, i, conj=c2, band=band, branch=b2) # right boundary position
		m = convert(PartialMPO, GTerm(pos1, pos2, coeff=1))
		A2 = m * cache.A
	
		k = GTEMPO.pos2pairindex(positions(m)[end])
		right = rightenv(cache, k)
		for tj in pos_j:k-1
			left = left * GrassmannTransferMatrix(tj, A2, cache.Bs...)
		end
		pos_j = k
		left′ = left * GrassmannTransferMatrix(k, A2, cache.Bs...)
		GFt[i] = GTEMPO.contract_center(left′, right) / Zvalue(cache)
	end

	return GFt
end