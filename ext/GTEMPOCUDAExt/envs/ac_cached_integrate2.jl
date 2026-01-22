
# set scaling for cache

cu_environments2(lattice::AbstractGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
						alg::IntegrationAlgorithm=ExactIntegrate()) = cuTwosideExpectationCache2(lattice, (A, B...))

function cuTwosideExpectationCache2(lattice::AbstractGrassmannLattice, As::Tuple; useHCache::Bool=DefaultUseCache)
	(ConjugationStyle(lattice.ordering) isa AdjacentConjugation) || throw(ArgumentError("only AdjacentConjugation supported for cached evaluation"))
	# (length(A) == length(B) == length(lattice)) || throw(DimensionMismatch())
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())
	Lhalf = div(length(lattice), 2)
	xs = As
	left = l_LL(xs...)
	hleft = useHCache ? CachedVector{typeof(left)}(undef, Lhalf+1) : Vector{typeof(left)}(undef, Lhalf+1)
	hleft_scaling = Vector{Float64}(undef, Lhalf+1)
	hleft_scaling[1] = _normalize!(left)
	hleft[1] = left
	left = tocu(left)
	for i in 1:Lhalf
		left = left * tocu(GrassmannTransferMatrix(i, xs...))
		hleft_scaling[i+1] = _normalize!(left)
		hleft[i+1] = fromcu(left)
		# hleft[i+1] = update_pair_left(hleft[i], i, xs..., trunc=trunc)
	end

	right = r_RR(xs...)
	hright = useHCache ? CachedVector{typeof(right)}(undef, Lhalf+1) : Vector{typeof(right)}(undef, Lhalf+1)
	hright_scaling = Vector{Float64}(undef, Lhalf+1)
	hright_scaling[Lhalf+1] = _normalize!(right)
	hright[Lhalf+1] = right
	right = tocu(right)
	for i in Lhalf:-1:1
		right = tocu(GrassmannTransferMatrix(i, xs...)) * right
		hright_scaling[i] = _normalize!(right)
		hright[i] = fromcu(right)
		# hright[i] = update_pair_right(hright[i+1], i, xs..., trunc=trunc)
	end
	return TwosideExpectationCache2(first(As), Base.tail(As), lattice, hleft, hleft_scaling, hright, hright_scaling)
end



# cached_gf_fast2.jl

# exhaust the first time index
function cu_cached_gf_fast(idx0::Int, lattice::AbstractGrassmannLattice, N::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; 
			cache::TwosideExpectationCache2=environments(lattice, A, Bs...),
			b1::Symbol, b2::Symbol, c1::Bool, c2::Bool, band::Int=1)
	pos1 = index(lattice, 2, conj=c1, band=band, branch=b1) 
	pos1′ = index(lattice, N, conj=c1, band=band, branch=b1) 
	pos2 = index(lattice, 1, conj=c2, band=band, branch=b2)
	if pos2 < min(pos1, pos1′)
		GFt = cu_cached_gf_fast_normal_order(idx0, lattice, N, A, Bs...; b1=b1, b2=b2, c1=c1, c2=c2, band=band, cache=cache)
	elseif pos2 > max(pos1, pos1′)
		GFt = cu_cached_gf_fast_reverse_order(idx0, lattice, N, A, Bs...; b1=b1, b2=b2, c1=c1, c2=c2, band=band, cache=cache)
	else
		error("time step indices should be monotonically ordered")
	end
	return GFt
end



function cu_cached_gf_fast_normal_order(idx0::Int, lattice::AbstractGrassmannLattice, idx1::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS};
					cache::TwosideExpectationCache2=environments(lattice, A, Bs...), 
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
	scaling::Float64 = 1.0
	for tj in j:k-1
		left = left * GrassmannTransferMatrix(tj, A2, cache.Bs...)
		scaling *= (_normalize!(left) / leftscaling(cache, tj+1))
	end
	pos_j = k
	left′ = left * GrassmannTransferMatrix(k, A2, cache.Bs...)
	scaling′ = scaling * (_normalize!(left′) / leftscaling(cache, k+1))
	GFt[2] = contract_center(left′, right) * scaling′ / contract_center(leftenv(cache, k+1), right)

	left = tocu(left)
	# calculated G1k
	for i in 3:N
		pos1 = index(lattice, idx0+i-1, conj=c1, band=band, branch=b1) # right boundary position
		pos2 = index(lattice, idx0+i-3, conj=c1, band=band, branch=b1)
		@assert pos2 < pos1
		m = convert(PartialMPO, GTerm(pos1, pos2, coeff=1))
		A2 = _mult_A(m, cache.A) 
	
		k = pos2pairindex(positions(m)[end])
		right = tocu(rightenv(cache, k))
		for tj in pos_j:k-1
			left = left * tocu(GrassmannTransferMatrix(tj, A2, cache.Bs...))
			scaling *= (_normalize!(left) / leftscaling(cache, tj+1))
		end
		pos_j = k
		left′ = left * tocu(GrassmannTransferMatrix(k, A2, cache.Bs...))
		scaling′ = scaling * (_normalize!(left′) / leftscaling(cache, k+1))
		GFt[i] = contract_center(left′, right) * scaling′ / contract_center(tocu(leftenv(cache, k+1)), right)
	end

	return GFt
end


function cu_cached_gf_fast_reverse_order(idx0::Int, lattice::AbstractGrassmannLattice, idx1::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS};
                    				 cache::TwosideExpectationCache2=environments(lattice, A, B...),
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
	scaling::Float64 = 1.0
	A2 = _mult_A(m, cache.A) 
	for tj in k:-1:j+1
		right = GrassmannTransferMatrix(tj, A2, cache.Bs...) * right
		scaling *= (_normalize!(right) / rightscaling(cache, tj-1))
	end
	pos_j = j
	right′ = GrassmannTransferMatrix(j, A2, cache.Bs...) * right 
	scaling′ = scaling * (_normalize!(right′) / rightscaling(cache, j-1))
	GFt[2] = contract_center(left, right′) * scaling′ / contract_center(left, rightenv(cache, j-1))

	right = tocu(right)
	for i in 3:N
		pos1 = index(lattice, idx0+i-1, conj=c1, band=band, branch=b1) # left boundary position
		pos2 = index(lattice, idx0+i-3, conj=c1, band=band, branch=b1)
		@assert pos1 < pos2
		m = convert(PartialMPO, GTerm(pos1, pos2, coeff=1))
		j = pos2pairindex(positions(m)[1])
		left = tocu(leftenv(cache, j))
		A2 = _mult_A(m, cache.A) 
		for tj in pos_j:-1:j+1
			right = tocu(GrassmannTransferMatrix(tj, A2, cache.Bs...)) * right
			scaling *= (_normalize!(right) / rightscaling(cache, tj-1))
		end
		pos_j = j
		right′ = tocu(GrassmannTransferMatrix(j, A2, cache.Bs...)) * right 
		scaling′ = scaling * (_normalize!(right′) / rightscaling(cache, j-1))
		GFt[i] = contract_center(left, right′) * scaling′ / contract_center(left, tocu(rightenv(cache, j-1)))
	
	end
	return GFt	
end


