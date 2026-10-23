# function contract_center(left::AbstractTensorMap{S, 1, 1}, right::AbstractTensorMap) where S 
# 	@tensor r = left[1,2] * right[2,1]
# 	return r
# end

function contract_center(left::AbstractParityTensorMap{<:Number, 1, N}, right::AbstractParityTensorMap{<:Number, N, 1}) where {N}
	@assert space(left, 1) == space(right, N+1)'

	s = space(left, 1) ← space(right, N+1)'
	r = fill!(similar(left, s), zero(scalartype(left)))

	cindA = ntuple(x->x+1, N)
	cindB = ntuple(x->N-x+1, N)
	contract!(r, left, ((1,), cindA), right, (cindB, (N+1,)), ((1,), (2,)), true, false, GrassmannBackend())
	return tr(r)
end


# DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS...) = error("l_LL not implemented for $(2+length(z)) mps")
# DMRG.r_RR(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS...) = error("r_RR not implemented for $(2+length(z)) mps")

# function DMRG.l_LL(x::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, one(vacuum), vacuum)
# 	return left
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum', vacuum)
# 	return left
# end

# function DMRG.r_RR(x::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum, one(vacuum))
# 	return left
# end

# function DMRG.r_RR(x::GrassmannMPS, y::GrassmannMPS)
# 	vacuum = space_l(x)
# 	right = isomorphism(Matrix{scalartype(x)}, vacuum, vacuum')
# 	return right
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum' ⊗ vacuum',  vacuum )
# 	return left
# end

# function DMRG.r_RR(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS)
# 	vacuum = space_l(x)
# 	right = isomorphism(Matrix{scalartype(x)}, vacuum ⊗ vacuum,  vacuum')
# 	return right
# end
# function contract_center(left::AbstractTensorMap{S, 2, 1}, right::AbstractTensorMap) where S 
# 	@tensor r = left[1,2,3] * right[3,2,1]
# 	return r
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS, u::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum' ⊗ vacuum' ,  vacuum ⊗ vacuum )
# 	return left
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS, u::GrassmannMPS, v::GrassmannMPS, w::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum' ⊗ vacuum' ⊗ vacuum',  vacuum ⊗ vacuum ⊗ vacuum )
# 	return left
# end

# function DMRG.l_LL(x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS, u::GrassmannMPS, v::GrassmannMPS, w::GrassmannMPS, r::GrassmannMPS)
# 	vacuum = space_l(x)
# 	left = isomorphism(Matrix{scalartype(x)}, vacuum' ⊗ vacuum' ⊗ vacuum' ⊗ vacuum',  vacuum ⊗ vacuum ⊗ vacuum )
# 	return left
# end

_mult_A(t, A::GrassmannMPS) = t * A
_mult_A(t, A::Vector{<:GrassmannMPS}) = [t * Aj for Aj in A]


# # the 2 and 3-th indices are physical ones to be traced out
# # nt = true means normal trace aa*
# function _trace_physical(a::AbstractTensorMap{S, 3, 1}; nt::Bool=true) where {S}
# 	r = TensorMap(ds->zeros(scalartype(a), ds), space(a, 1) ← space(a, 4)' )
# 	I = sectortype(a)
# 	for (f1, f2) in fusiontrees(a)
# 		if f1.uncoupled[2] == f1.uncoupled[3]
# 			f0 = FusionTree{I}((f1.uncoupled[1],), f1.coupled, (f1.isdual[1],), ())
# 			coef = (isodd(f1.uncoupled[2].n) && (!nt)) ? -1 : 1
# 			@tensor r[f0, f2][1,3] += coef * a[f1, f2][1,2,2,3]
# 			# axpy!(1, tmp, r[f0, f2]) 
# 		end
# 	end
# 	return r
# end


# function swap12!(twositemps::AbstractTensorMap{S, 3, 1}) where S
# 	for (f1, f2) in fusiontrees(twositemps)
# 		if isodd(f1.uncoupled[1].n) && isodd(f1.uncoupled[2].n)
# 			twositemps[f1, f2] .*= -1
# 		end
# 	end	
# 	return twositemps
# end


# function _fuse_physical(m::AbstractTensorMap{S, 3, 1}) where S
# 	cod = ProductSpace{S}((space(m, 1), space(m, 2)))
# 	dom = ProductSpace{S}((space(m, 4)',))
# 	r = TensorMap(ds->zeros(scalartype(m), ds), cod ← dom)
# 	for (f1, f2) in fusiontrees(m)
# 		if f1.uncoupled[2] == f1.uncoupled[3]
# 			if iseven(f1.uncoupled[2].n)
# 				f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2]), f1.coupled, (f1.isdual[1],f1.isdual[2]))
# 				r[f1′, f2] .+= m[f1, f2][:, :, 1, :]
# 			end
# 		else
# 			isdual2 = isodd(f1.uncoupled[2].n) ? f1.isdual[2] : f1.isdual[3]
# 			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1)), f1.coupled, (f1.isdual[1], isdual2))
# 			r[f1′, f2] .+= m[f1, f2][:, :, 1, :]
# 		end
# 	end
# 	return r
# end


# function swap34!(twositemps::AbstractTensorMap{S, 1, 3}) where S
# 	for (f1, f2) in fusiontrees(twositemps)
# 		if isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[3].n)
# 			twositemps[f1, f2] .*= -1
# 		end
# 	end	
# 	return twositemps
# end

# function _fuse_physical(m::AbstractTensorMap{S, 1, 3}) where S
# 	cod = ProductSpace{S}((space(m, 1),))
# 	dom = ProductSpace{S}((space(m, 3)', space(m, 4)'))
# 	r = TensorMap(ds->zeros(scalartype(m), ds), cod ← dom)
# 	for (f1, f2) in fusiontrees(m)
# 		if f2.uncoupled[1] == f2.uncoupled[2]
# 			if iseven(f2.uncoupled[1].n)
# 				f2′ = FusionTree((f2.uncoupled[2],f2.uncoupled[3]), f2.coupled, (f2.isdual[2],f2.isdual[3]))
# 				r[f1, f2′] .+= m[f1, f2][:, 1, :, :]
# 			end
# 		else
# 			isdual2 = isodd(f2.uncoupled[2].n) ? f2.isdual[2] : f2.isdual[1]
# 			f2′ = FusionTree((Z2Irrep(1), f2.uncoupled[3]), f2.coupled, (isdual2, f2.isdual[3]))
# 			r[f1, f2′] .+= m[f1, f2][:, 1, :, :]
# 		end
# 	end
# 	return r
# end

# cached integration
abstract type AbstractExpectationCache end

"""
	Zvalue(x::AbstractExpectationCache)

Return the value of partition function
"""
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

"""
    environments(lattice::AbstractGrassmannLattice, A, B...; alg=ExactIntegrate())

Precompute the left/right environment tensors for the product `A * B * ...` on
`lattice` and return them as an `AbstractExpectationCache`. The cache carries
the partition function (`Zvalue`) and allows cheap evaluation of local
observables and Green's functions through `expectationvalue` / `cached_gf`
without repeating the contraction.
"""
environments(lattice::AbstractGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
						alg::IntegrationAlgorithm=ExactIntegrate()) = _environments(alg, lattice, A, B...)

# cached_integrate_util(lattice::AbstractGrassmannLattice, j::Int, k::Int, cache::AbstractExpectationCache, A::Union{GrassmannMPS, Vector}, Bs::GrassmannMPS...; kwargs...) = _cached_integrate_util(
# 	lattice, j, k, cache, A, Bs...; kwargs...) / Zvalue(cache)

expectationvalue(m::Union{GTerm, ExpGTerm}, cache::AbstractExpectationCache) = expectationvalue(convert(PartialMPO, m), cache)

struct TwosideExpectationCache{M<:GrassmannMPS, G<:Tuple, L<:AbstractGrassmannLattice, Hl, Hr} <: AbstractExpectationCache
	A::M
	Bs::G
	lattice::L
	hleft::Hl
	hright::Hr	
end
rightenv(x::TwosideExpectationCache, j::Int) = x.hright[j+1]
Base.length(x::TwosideExpectationCache) = length(x.lattice)

function TwosideExpectationCache(lattice::AbstractGrassmannLattice, As::Tuple)
	(ConjugationStyle(lattice.ordering) isa AdjacentConjugation) || throw(ArgumentError("only AdjacentConjugation supported for cached evaluation"))
	# (length(A) == length(B) == length(lattice)) || throw(DimensionMismatch())
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())
	Lhalf = div(length(lattice), 2)
	xs = As
	left = l_LL(xs...)
	hleft = Vector{typeof(left)}(undef, Lhalf+1)
	hleft[1] = left
	for i in 1:Lhalf
		left = left * GrassmannTransferMatrix(i, xs...)
		hleft[i+1] = left
		# hleft[i+1] = update_pair_left(hleft[i], i, xs..., trunc=trunc)
	end

	right = r_RR(xs...)
	hright = Vector{typeof(right)}(undef, Lhalf+1)
	hright[Lhalf+1] = right
	for i in Lhalf:-1:1
		right = GrassmannTransferMatrix(i, xs...) * right
		hright[i] = right
		# hright[i] = update_pair_right(hright[i+1], i, xs..., trunc=trunc)
	end
	return TwosideExpectationCache(first(As), Base.tail(As), lattice, hleft, hright)
end

_environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}) = TwosideExpectationCache(lattice, (A, B...))
# _environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::GrassmannMPS, B::GrassmannMPS, C::GrassmannMPS; kwargs...) = TwosideExpectationCache(lattice, (A, B, C); kwargs...)


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
expectationvalue(m::PartialMPO, cache::TwosideExpectationCache) = expectation(m, cache) / Zvalue(cache)
function expectation(m::PartialMPO, cache::TwosideExpectationCache)
	j, k = positions(m)[1], positions(m)[end]
	j, k = pos2pairindex(j), pos2pairindex(k)
	left = leftenv(cache, j)  
	right = rightenv(cache, k) 
	A2 = _mult_A(m, cache.A)
	for tj in k:-1:j
		right = GrassmannTransferMatrix(tj, A2, cache.Bs...) * right 
		# right = update_pair_right(right, tj, A2, cache.Bs..., trunc=trunc)
	end	
	return contract_center(left, right) 	
end

# 4 mps
struct LeftExpectationCache{M<:GrassmannMPS, G<:Tuple, L<:AbstractGrassmannLattice, H} <: AbstractExpectationCache
	A::M
	Bs::G
	lattice::L
	hleft::H
end
Base.length(x::LeftExpectationCache) = length(x.lattice)
function LeftExpectationCache(lattice::AbstractGrassmannLattice, As::Tuple)
	(ConjugationStyle(lattice.ordering) isa AdjacentConjugation) || throw(ArgumentError("only AdjacentConjugation supported for cached evaluation"))
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())
	Lhalf = div(length(lattice), 2)
	xs = As
	left = l_LL(xs...)
	hleft = Vector{typeof(left)}(undef, Lhalf+1)
	hleft[1] = left
	for i in 1:Lhalf
		hleft[i+1] = hleft[i] * GrassmannTransferMatrix(i, xs...)
		# hleft[i+1] = update_pair_left(hleft[i], i, xs..., trunc=trunc)
	end
	return LeftExpectationCache(first(As), Base.tail(As), lattice, hleft)
end
# _environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::GrassmannMPS, B::GrassmannMPS, C::GrassmannMPS, D::GrassmannMPS, E::GrassmannMPS...; 
# 						kwargs...) = LeftExpectationCache(lattice, (A, B, C, D, E...); kwargs...)

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
expectationvalue(m::PartialMPO, cache::LeftExpectationCache) = expectation(m, cache) / Zvalue(cache)
function expectation(m::PartialMPO, cache::LeftExpectationCache)
	j, k = positions(m)[1], positions(m)[end]
	j = pos2pairindex(j)
	left = leftenv(cache, j)  
	A2 = _mult_A(m, cache.A)
	for tj in j:pos2pairindex(length(cache))
		left = left * GrassmannTransferMatrix(tj, A2, cache.Bs...)
		# left = update_pair_left(left, tj, A2, cache.Bs..., trunc=trunc)
	end	
	return TK.scalar(left) 
end

struct VectorExpectationCache{C<:AbstractExpectationCache} <: AbstractExpectationCache
	caches::Vector{C}
end
Base.length(x::VectorExpectationCache) = length(x.caches[1])
Zvalue(x::VectorExpectationCache) = sum(Zvalue.(x.caches))

function _environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::Vector{<:GrassmannMPS}, B::Vararg{GrassmannMPS})
	caches = [_environments(alg, lattice, Aj, B...) for Aj in A]
	return VectorExpectationCache(caches)
end

# function _cached_integrate_util(lattice::AbstractGrassmannLattice, j::Int, k::Int, cache::VectorExpectationCache, A2::Vector{<:GrassmannMPS}, Bs::GrassmannMPS...; kwargs...)
# 	@assert length(cache.caches) == length(A2)
# 	return sum([_cached_integrate_util(lattice, j, k, cj, Aj, Bs...; kwargs...) for (cj, Aj) in zip(cache.caches, A2)])
# end
expectation(m::PartialMPO, cache::VectorExpectationCache) = sum(expectation(m, cj) for cj in cache.caches)
expectationvalue(m::PartialMPO, cache::VectorExpectationCache) = expectation(m, cache) / Zvalue(cache)

