cached_Gτ_fast(lattice::ImagGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; c1::Bool=false, c2::Bool=true, kwargs...) = cached_gf_fast(
				lattice, A, B...;  c1=c1, c2=c2, b1=:τ, b2=:τ, kwargs...)
function cached_Gτ_fast(lattice::MixedGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; c1::Bool=false, c2::Bool=true, kwargs...)
	r = cached_gf_fast(lattice, A, B...;  c1=c1, c2=c2, b1=:τ, b2=:τ, kwargs...)
	r[end] = 1 - r[1]
	return r
end 

cached_Gt_fast(lattice::RealGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; c1::Bool=true, c2::Bool=false, kwargs...) = cached_gf_fast(
				lattice, A, B...; c1=c1, c2=c2, kwargs...)
cached_Gm_fast(lattice::MixedGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}; c1::Bool=true, c2::Bool=false, kwargs...) = cached_gf_fast(
				lattice, A, B...; c1=c1, c2=c2, kwargs...)

cached_greater_fast(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, A::GrassmannMPS, B::Vararg{GrassmannMPS}; kwargs...) = cached_gf_fast(
	lattice, A, B...; b1=:+, b2=:+, c1=false, c2=true, kwargs...)
cached_lesser_fast(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, A::GrassmannMPS, B::Vararg{GrassmannMPS}; kwargs...) = -cached_gf_fast(
	lattice, A, B...; b1=:+, b2=:-, c1=false, c2=true, kwargs...)



"""
	cached_gf_fast(lattice::Union{ImagGrassmannLattice, RealGrassmannLattice}, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; kwargs...)

Similar to cached_gf, but use a even faster version when calculating 
many ⟨aᵢa₁⟩s with differen distances |i-1|, by more cleverly using the cache

Currently, this function may contain errors if b1=:- and b2=:+, so one should 
avoid this situation
In the future, we may either solve this problem or simply 
throw an error when this situation is met
"""
function cached_gf_fast(lattice::Union{ImagGrassmannLattice, RealGrassmannLattice}, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; kwargs...)
	GFt = cached_gf_fast(lattice, lattice.k, A, Bs...; kwargs...)
	if lattice isa ImagGrassmannLattice
		GFt[end] = 1 - GFt[1]
	end
	return GFt
end 
function cached_gf_fast(lattice::MixedGrassmannLattice, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; b1::Symbol, b2::Symbol,kwargs...)
	if (b1 == :+) || (b1 == :-)
		return cached_gf_fast(lattice, lattice.kt, A, Bs...; b1=b1, b2=b2, kwargs...)
	else
		return cached_gf_fast(lattice, lattice.kτ, A, Bs...; b1=b1, b2=b2, kwargs...)
	end
end
function cached_gf_fast(lattice::AbstractGrassmannLattice, N::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; b1::Symbol, b2::Symbol, c1::Bool, c2::Bool, band::Int=1, kwargs...)
	pos1 = index(lattice, 2, conj=c1, band=band, branch=b1) 
	pos1′ = index(lattice, N, conj=c1, band=band, branch=b1) 
	pos2 = index(lattice, 1, conj=c2, band=band, branch=b2)
	if pos2 < min(pos1, pos1′)
		GFt = cached_gf_fast_normal_order(lattice, N, A, Bs...; b1=b1, b2=b2, c1=c1, c2=c2, band=band, kwargs...)
	elseif pos2 > max(pos1, pos1′)
		GFt = cached_gf_fast_reverse_order(lattice, N, A, Bs...; b1=b1, b2=b2, c1=c1, c2=c2, band=band, kwargs...)
	else
		error("time step indices should be monotonically ordered")
	end
	return GFt
end

# function cached_gf_fast(lattice::Union{ImagGrassmannLattice, RealGrassmannLattice}, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; kwargs...)
# 	if TimeOrderingStyle(lattice) isa TimeAscending
# 		GFt = cached_gf_fast_normal_order(lattice, lattice.k, A, Bs...; kwargs...)
# 	else
# 		GFt = cached_gf_fast_reverse_order(lattice, lattice.k, A, Bs...; kwargs...)
# 	end
# 	if lattice isa ImagGrassmannLattice
# 		GFt[end] = 1 - GFt[1]
# 	end
# 	return GFt
# end

# function cached_gf_fast(lattice::MixedGrassmannLattice, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; b1::Symbol, b2::Symbol,kwargs...)
# 	if (b2 == :+) || (b2 == :-)
# 		if RealTimeOrderingStyle(lattice) isa TimeAscending
# 			GFt = cached_gf_fast_normal_order(lattice, lattice.kt, A, Bs...; b1=b1, b2=b2, kwargs...)
# 		else
# 			GFt = cached_gf_fast_reverse_order(lattice, lattice.kt, A, Bs...; b1=b1, b2=b2, kwargs...)
# 		end		
# 		return GFt
# 	else
# 		if ImaginaryTimeOrderingStyle(lattice) isa TimeAscending
# 			GFt = cached_gf_fast_normal_order(lattice, lattice.kτ, A, Bs...; b1=b1, b2=b2, kwargs...)
# 		else
# 			GFt = cached_gf_fast_reverse_order(lattice, lattice.kτ, A, Bs...; b1=b1, b2=b2, kwargs...)
# 		end		
# 		return GFt
# 	end
# end

function cached_gf_fast_normal_order(lattice::AbstractGrassmannLattice, N::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS}; 
					cache::TwosideExpectationCache=environments(lattice, A, Bs...), 
					b1::Symbol, b2::Symbol, c1::Bool, c2::Bool, band::Int=1, kwargs...)
	@assert N > 0
	# calculated G11
	a, b = ContourIndex(1, conj=c1, branch=b1, band=band), ContourIndex(1, conj=c2, branch=b2, band=band)
	g0 = cached_gf(lattice, a, b, A, Bs...; cache=cache, kwargs...)
	GFt = Vector{scalartype(A)}(undef, N)
	GFt[1] = g0
	(N == 1) && return GFt
	
	# calculated G12
	pos2 = index(lattice, 1, conj=c2, band=band, branch=b2) # left boundary position
	pos1 = index(lattice, 2, conj=c1, band=band, branch=b1) # right boundary position
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
		pos1 = index(lattice, i, conj=c1, band=band, branch=b1) # right boundary position
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


function cached_gf_fast_reverse_order(lattice::AbstractGrassmannLattice, N::Int, A::GrassmannMPS, Bs::Vararg{GrassmannMPS};
                    				 cache::TwosideExpectationCache=environments(lattice, A, B...),
                    				 b1::Symbol, b2::Symbol, c1::Bool, c2::Bool, band::Int=1, kwargs...)
	@assert N > 0
	a, b = ContourIndex(1, conj=c1, branch=b1, band=band), ContourIndex(1, conj=c2, branch=b2, band=band)
	g0 = cached_gf(lattice, a, b, A, Bs...; cache=cache, kwargs...)
	GFt = Vector{scalartype(A)}(undef, N)
	GFt[1] = g0
	(N == 1) && return GFt

	pos1 = index(lattice, 2, conj=c1, band=band, branch=b1) # left boundary position
	pos2 = index(lattice, 1, conj=c2, band=band, branch=b2) # right boundary position
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
		pos1 = index(lattice, i, conj=c1, band=band, branch=b1) # left boundary position
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