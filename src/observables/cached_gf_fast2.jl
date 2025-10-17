
# add idx0 for RealGrassmannLattice. based on cached_gf_fast.jl


cached_Gt_fast(idx0::Int, lattice::RealGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; c1::Bool=true, c2::Bool=false, kwargs...) = cached_gf_fast(
				idx0, lattice, A, B...; c1=c1, c2=c2, kwargs...)

cached_greater_fast(idx0::Int, lattice::RealGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; kwargs...) = cached_gf_fast(
	idx0, lattice, A, B...; b1=:+, b2=:+, c1=false, c2=true, kwargs...)
cached_lesser_fast(idx0::Int, lattice::RealGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; kwargs...) = -cached_gf_fast(
	idx0, lattice, A, B...; b1=:+, b2=:-, c1=false, c2=true, kwargs...)

cached_greater_fast(idx0::Int, lattice::RealGrassmannLattice, N::Int, A::GrassmannMPS, B::Vararg{GrassmannMPS}; kwargs...) = cached_gf_fast(
	idx0, lattice, N, A, B...; b1=:+, b2=:+, c1=false, c2=true, kwargs...)
cached_lesser_fast(idx0::Int, lattice::RealGrassmannLattice, N::Int, A::GrassmannMPS, B::Vararg{GrassmannMPS}; kwargs...) = -cached_gf_fast(
	idx0, lattice, N, A, B...; b1=:+, b2=:-, c1=false, c2=true, kwargs...)


function cached_gf_fast(idx0::Int, lattice::RealGrassmannLattice, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; kwargs...)
	cached_gf_fast(idx0, lattice, lattice.k, A, Bs...; kwargs...)
end 

# exhaust the first time index
function cached_gf_fast(idx0::Int, lattice::AbstractGrassmannLattice, N::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; 
			cache::TwosideExpectationCache=environments(lattice, A, Bs...),
			b1::Symbol, b2::Symbol, c1::Bool, c2::Bool, band::Int=1)
	pos1 = index(lattice, 2, conj=c1, band=band, branch=b1) 
	pos1′ = index(lattice, N, conj=c1, band=band, branch=b1) 
	pos2 = index(lattice, 1, conj=c2, band=band, branch=b2)
	if pos2 < min(pos1, pos1′)
		GFt = cached_gf_fast_normal_order(idx0, lattice, N, A, Bs...; b1=b1, b2=b2, c1=c1, c2=c2, band=band, cache=cache)
	elseif pos2 > max(pos1, pos1′)
		GFt = cached_gf_fast_reverse_order(idx0, lattice, N, A, Bs...; b1=b1, b2=b2, c1=c1, c2=c2, band=band, cache=cache)
	else
		error("time step indices should be monotonically ordered")
	end
	return GFt
end



function cached_gf_fast_normal_order(idx0::Int, lattice::AbstractGrassmannLattice, idx1::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS};
					cache::TwosideExpectationCache=environments(lattice, A, Bs...), 
					b1::Symbol, b2::Symbol, c1::Bool, c2::Bool, band::Int=1, kwargs...)
	@assert idx1 > idx0

	N = idx1 - idx0 + 1
	# calculated G11
	a, b = ContourIndex(idx0, conj=c1, branch=b1, band=band), ContourIndex(idx0, conj=c2, branch=b2, band=band)
	g0 = cached_gf(lattice, (a, b), A, Bs...; cache=cache, kwargs...)
	GFt = Vector{scalartype(A)}(undef, N)
	GFt[1] = g0
	(N == 1) && return GFt
	
	# calculated G12
	pos2 = index(lattice, idx0, conj=c2, band=band, branch=b2) # left boundary position
	pos1 = index(lattice, idx0+1, conj=c1, band=band, branch=b1) # right boundary position
	@assert pos2 < pos1
	m = convert(PartialMPO, GTerm(pos1, pos2, coeff=1))
	A2 = _mult_A(m, cache.A) 

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
	for i in 3:N
		pos1 = index(lattice, idx0+i-1, conj=c1, band=band, branch=b1) # right boundary position
		pos2 = index(lattice, idx0+i-3, conj=c1, band=band, branch=b1)
		@assert pos2 < pos1
		m = convert(PartialMPO, GTerm(pos1, pos2, coeff=1))
		A2 = _mult_A(m, cache.A) 
	
		k = pos2pairindex(positions(m)[end])
		right = rightenv(cache, k)
		for tj in pos_j:k-1
			left = left * GrassmannTransferMatrix(tj, A2, cache.Bs...)
		end
		pos_j = k
		left′ = left * GrassmannTransferMatrix(k, A2, cache.Bs...)
		GFt[i] = contract_center(left′, right) / Zvalue(cache)
	end

	return GFt
end


function cached_gf_fast_reverse_order(idx0::Int, lattice::AbstractGrassmannLattice, idx1::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS};
                    				 cache::TwosideExpectationCache=environments(lattice, A, B...),
                    				 b1::Symbol, b2::Symbol, c1::Bool, c2::Bool, band::Int=1, kwargs...)
	@assert idx1 > idx0

	N = idx1 - idx0 + 1
	a, b = ContourIndex(idx0, conj=c1, branch=b1, band=band), ContourIndex(idx0, conj=c2, branch=b2, band=band)
	g0 = cached_gf(lattice, (a, b), A, Bs...; cache=cache, kwargs...)
	GFt = Vector{scalartype(A)}(undef, N)
	GFt[1] = g0
	(N == 1) && return GFt

	pos1 = index(lattice, idx0+1, conj=c1, band=band, branch=b1) # left boundary position
	pos2 = index(lattice, idx0, conj=c2, band=band, branch=b2) # right boundary position
	@assert pos1 < pos2
	m = convert(PartialMPO, GTerm(pos1, pos2, coeff=1))
	j, k = pos2pairindex(positions(m)[1]), pos2pairindex(positions(m)[end])
	left = leftenv(cache, j)
	right = rightenv(cache, k)
	A2 = _mult_A(m, cache.A) 
	for tj in k:-1:j+1
		right = GrassmannTransferMatrix(tj, A2, cache.Bs...) * right
	end
	pos_j = j
	right′ = GrassmannTransferMatrix(j, A2, cache.Bs...) * right 
	GFt[2] = contract_center(left, right′) / Zvalue(cache)

	for i in 3:N
		pos1 = index(lattice, idx0+i-1, conj=c1, band=band, branch=b1) # left boundary position
		pos2 = index(lattice, idx0+i-3, conj=c1, band=band, branch=b1)
		@assert pos1 < pos2
		m = convert(PartialMPO, GTerm(pos1, pos2, coeff=1))
		j = pos2pairindex(positions(m)[1])
		left = leftenv(cache, j)
		A2 = _mult_A(m, cache.A) 
		for tj in pos_j:-1:j+1
			right = GrassmannTransferMatrix(tj, A2, cache.Bs...) * right
		end
		pos_j = j
		right′ = GrassmannTransferMatrix(j, A2, cache.Bs...) * right 
		GFt[i] = contract_center(left, right′) / Zvalue(cache)
	
	end
	return GFt	
end

