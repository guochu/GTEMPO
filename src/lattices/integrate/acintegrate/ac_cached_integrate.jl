# cached integration
abstract type AbstractExpectationCache end

Zvalue(x::AbstractExpectationCache) = TK.scalar(x.hleft[end])
"""
	leftenv(x::AbstractExpectationCache, j::Int)

Cache at the left of the j-th time step
"""
leftenv(x::AbstractExpectationCache, j::Int) = x.hleft[j]
"""
	rightenv(x::AbstractExpectationCache, j::Int)

Cache at the right of the j-th time step
"""
rightenv(x::AbstractExpectationCache, j::Int) = error("rightenv not implemented for cache type $(typeof(x))")

DMRG.environments(lattice::AbstractGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS, C::GrassmannMPS...; 
						alg::IntegrationAlgorithm=ExactIntegrate(), kwargs...) = _environments(alg, lattice, A, B, C...; kwargs...)

# cached_integrate_util(lattice::AbstractGrassmannLattice, j::Int, k::Int, cache::AbstractExpectationCache, A::Union{GrassmannMPS, Vector}, Bs::GrassmannMPS...; kwargs...) = _cached_integrate_util(
# 	lattice, j, k, cache, A, Bs...; kwargs...) / Zvalue(cache)

integrate(m::Union{GTerm, ExpGTerm}, cache::AbstractExpectationCache; kwargs...) = integrate(convert(PartialMPO, m), cache; kwargs...)

struct TwosideExpectationCache{M<:GrassmannMPS, G<:Tuple, L<:AbstractGrassmannLattice, Hl, Hr} <: AbstractExpectationCache
	A::M
	Bs::G
	lattice::L
	hleft::Hl
	hright::Hr	
end
rightenv(x::TwosideExpectationCache, j::Int) = x.hright[j+1]
Base.length(x::TwosideExpectationCache) = length(x.lattice)

function TwosideExpectationCache(lattice::AbstractGrassmannLattice, As::Tuple; trunc::TruncationScheme=DefaultIntegrationTruncation) 
	(ConjugationStyle(lattice.ordering) isa AdjacentConjugation) || throw(ArgumentError("only AdjacentConjugation supported for cached evaluation"))
	# (length(A) == length(B) == length(lattice)) || throw(DimensionMismatch())
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())
	Lhalf = div(length(lattice), 2)
	xs = As
	left = _l_LL(xs...)
	hleft = Vector{typeof(left)}(undef, Lhalf+1)
	hleft[1] = left
	for i in 1:Lhalf
		hleft[i+1] = update_pair_left(hleft[i], i, xs..., trunc=trunc)
	end

	right = _r_RR(xs...)
	hright = Vector{typeof(right)}(undef, Lhalf+1)
	hright[Lhalf+1] = right
	for i in Lhalf:-1:1
		hright[i] = update_pair_right(hright[i+1], i, xs..., trunc=trunc)
	end
	return TwosideExpectationCache(first(As), Base.tail(As), lattice, hleft, hright)
end

_environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::GrassmannMPS, B::GrassmannMPS; kwargs...) = TwosideExpectationCache(lattice, (A, B); kwargs...)
_environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::GrassmannMPS, B::GrassmannMPS, C::GrassmannMPS; kwargs...) = TwosideExpectationCache(lattice, (A, B, C); kwargs...)


# function _cached_integrate_util(lattice::AbstractGrassmannLattice, j::Int, k::Int, cache::TwosideExpectationCache, A2::GrassmannMPS, Bs::GrassmannMPS...; 
# 								trunc::TruncationScheme=DefaultIntegrationTruncation)
# 	@assert cache.lattice === lattice
# 	@assert all(v -> v[1] === v[2], zip(cache.Bs, Bs)) || throw(DimensionMismatch())
# 	@assert j <= k
# 	j, k = pos2pairindex(j), pos2pairindex(k)
# 	B2 = Bs
# 	left = leftenv(cache, j)  
# 	right = rightenv(cache, k) 
# 	for tj in k:-1:j
# 		right = update_pair_right(right, tj, A2, B2..., trunc=trunc)
# 	end
# 	r = contract_center(left, right) 
# 	return r 
# end
integrate(m::PartialMPO, cache::TwosideExpectationCache; trunc::TruncationScheme=DefaultIntegrationTruncation) = integrate_unnormalized(
			m, cache, trunc=trunc) / Zvalue(cache)
function integrate_unnormalized(m::PartialMPO, cache::TwosideExpectationCache; trunc::TruncationScheme=DefaultIntegrationTruncation)
	j, k = positions(m)[1], positions(m)[end]
	j, k = pos2pairindex(j), pos2pairindex(k)
	left = leftenv(cache, j)  
	right = rightenv(cache, k) 
	A2 = _mult_A(m, cache.A)
	for tj in k:-1:j
		right = update_pair_right(right, tj, A2, cache.Bs..., trunc=trunc)
	end	
	return contract_center(left, right) 	
end
# cached_integrate_util(lattice::AbstractGrassmannLattice, j::Int, k::Int, cache::TwosideExpectationCache, A2::GrassmannMPS, Bs::GrassmannMPS...; kwargs...) = _cached_integrate_util(
# 	lattice, j, k, cache, A2, Bs...; kwargs...) / Zvalue(cache)

# 4 mps
struct LeftExpectationCache{M<:GrassmannMPS, G<:Tuple, L<:AbstractGrassmannLattice, H} <: AbstractExpectationCache
	A::M
	Bs::G
	lattice::L
	hleft::H
end
Base.length(x::LeftExpectationCache) = length(x.lattice)
function LeftExpectationCache(lattice::AbstractGrassmannLattice, As::Tuple; trunc::TruncationScheme=DefaultIntegrationTruncation)
	(ConjugationStyle(lattice.ordering) isa AdjacentConjugation) || throw(ArgumentError("only AdjacentConjugation supported for cached evaluation"))
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())
	Lhalf = div(length(lattice), 2)
	xs = As
	left = _l_LL(xs...)
	hleft = Vector{typeof(left)}(undef, Lhalf+1)
	hleft[1] = left
	for i in 1:Lhalf
		hleft[i+1] = update_pair_left(hleft[i], i, xs..., trunc=trunc)
	end
	return LeftExpectationCache(first(As), Base.tail(As), lattice, hleft)
end
_environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::GrassmannMPS, B::GrassmannMPS, C::GrassmannMPS, D::GrassmannMPS, E::GrassmannMPS...; 
						kwargs...) = LeftExpectationCache(lattice, (A, B, C, D, E...); kwargs...)

# A2 is temporary and overwritten
# function _cached_integrate_util(lattice::AbstractGrassmannLattice, j::Int, k::Int, cache::LeftExpectationCache, A2::GrassmannMPS, Bs::GrassmannMPS...; 
# 								trunc::TruncationScheme=DefaultIntegrationTruncation)
# 	# @assert (cache.B === B) && (cache.C === C) && (cache.D === D) && (cache.lattice === lattice) 
# 	@assert cache.lattice === lattice
# 	@assert all(v -> v[1] === v[2], zip(cache.Bs, Bs)) || throw(DimensionMismatch())
# 	@assert j <= k
# 	j = pos2pairindex(j)
# 	B2 = Bs
# 	left = leftenv(cache, j)  
# 	for tj in j:pos2pairindex(length(lattice))
# 		left = update_pair_left(left, tj, A2, B2..., trunc=trunc)
# 	end
# 	return TK.scalar(left) 
# end
# cached_integrate_util(lattice::AbstractGrassmannLattice, j::Int, k::Int, cache::LeftExpectationCache, A2::GrassmannMPS, Bs::GrassmannMPS...; kwargs...) = _cached_integrate_util(
# 	lattice, j, k, cache, A2, Bs...; kwargs...) / Zvalue(cache)
integrate(m::PartialMPO, cache::LeftExpectationCache; trunc::TruncationScheme=DefaultIntegrationTruncation) = integrate_unnormalized(
			m, cache, trunc=trunc) / Zvalue(cache)
function integrate_unnormalized(m::PartialMPO, cache::LeftExpectationCache; trunc::TruncationScheme=DefaultIntegrationTruncation)
	j, k = positions(m)[1], positions(m)[end]
	j = pos2pairindex(j)
	left = leftenv(cache, j)  
	A2 = _mult_A(m, cache.A)
	for tj in j:pos2pairindex(length(cache))
		left = update_pair_left(left, tj, A2, cache.Bs..., trunc=trunc)
	end	
	return TK.scalar(left) 
end

struct VectorExpectationCache{C<:AbstractExpectationCache} <: AbstractExpectationCache
	caches::Vector{C}
end
Base.length(x::VectorExpectationCache) = length(x.caches[1])
Zvalue(x::VectorExpectationCache) = sum(Zvalue.(x.caches))

function _environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; kwargs...)
	caches = [_environments(alg, lattice, Aj, B...; kwargs...) for Aj in A]
	return VectorExpectationCache(caches)
end

# function _cached_integrate_util(lattice::AbstractGrassmannLattice, j::Int, k::Int, cache::VectorExpectationCache, A2::Vector{<:GrassmannMPS}, Bs::GrassmannMPS...; kwargs...)
# 	@assert length(cache.caches) == length(A2)
# 	return sum([_cached_integrate_util(lattice, j, k, cj, Aj, Bs...; kwargs...) for (cj, Aj) in zip(cache.caches, A2)])
# end

function integrate(m::PartialMPO, cache::VectorExpectationCache; trunc::TruncationScheme=DefaultIntegrationTruncation)
	return sum(integrate_unnormalized(m, cj; trunc=trunc) for cj in cache.caches) / Zvalue(cache)
end


