
function cached_nn2(lattice::AbstractGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol=:τ, b2::Symbol=b1, 
                    band::Union{Int, Tuple{Int, Int}}=1)
	if isa(band, Int)
		band = (band, band)
	end
	band1, band2 = band
	a1, a2 = get_nn_contour_pos(lattice, i, band1, b1)
	if (i == j) && (b1 == b2) && (band1==band2)
		return cached_gf(lattice, (a1, a2), A, B...; cache=cache)
	end
	a3, a4 = get_nn_contour_pos(lattice, j, band2, b2)
	return cached_gf(lattice, (a1, a2, a3, a4), A, B...; cache=cache)
end

function cached_nn(lattice::AbstractGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol=:τ, b2::Symbol=b1, 
                    band::Union{Int, Tuple{Int, Int}}=1)
	if isa(band, Int)
		band = (band, band)
	end
	band1, band2 = band
	a, b = ContourIndex(i, conj=false, band=band1, branch=b1), ContourIndex(j, conj=false, band=band2, branch=b2)
	return _nn(a, b, cache)
end

function _nn(a::ContourIndex, b::ContourIndex, cache::TwosideExpectationCache)
	j, k = cache.lattice[a], cache.lattice[b]
	j, k = pos2pairindex(j), pos2pairindex(k)
	if j > k
		j, k = k, j
		a, b = b, a
	end
	left = leftenv(cache, j)  
	right = rightenv(cache, k) 
	if a == b
		A2 = _insert_n(cache.lattice, cache.A, a.j, branch(a), a.band)
	else
		A2 = _insert_nn(cache.lattice, cache.A, a.j, b.j, branch(a), branch(b), a.band, b.band)
	end
	for tj in k:-1:j
		right = GrassmannTransferMatrix(tj, A2, cache.Bs...) * right 
		# right = update_pair_right(right, tj, A2, cache.Bs..., trunc=trunc)
	end	
	return contract_center(left, right) / Zvalue(cache)	
end