
cu_environments(lattice::AbstractGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
						alg::IntegrationAlgorithm=ExactIntegrate()) = _cu_environments(alg, lattice, A, B...)

function cuTwosideExpectationCache(lattice::AbstractGrassmannLattice, As::Tuple; useHCache::Bool=DefaultUseCache) 
	(ConjugationStyle(lattice.ordering) isa AdjacentConjugation) || throw(ArgumentError("only AdjacentConjugation supported for cached evaluation"))
	# (length(A) == length(B) == length(lattice)) || throw(DimensionMismatch())
	(all(v->length(v)==length(lattice), As)) || throw(DimensionMismatch())
	Lhalf = div(length(lattice), 2)
	xs = As
	left = l_LL(xs...)
	hleft = useHCache ? CachedVector{typeof(left)}(undef, Lhalf+1) : Vector{typeof(left)}(undef, Lhalf+1)
	hleft[1] = left
	left = tocu(left)
	for i in 1:Lhalf
		left = left * tocu(GrassmannTransferMatrix(i, xs...))
		hleft[i+1] = fromcu(left)
	end

	right = r_RR(xs...)
	hright = useHCache ? CachedVector{typeof(right)}(undef, Lhalf+1) : Vector{typeof(right)}(undef, Lhalf+1)
	hright[Lhalf+1] = right
	right = tocu(right)
	for i in Lhalf:-1:1
		right = tocu(GrassmannTransferMatrix(i, xs...)) * right
		hright[i] = fromcu(right)
	end
	return TwosideExpectationCache(first(As), Base.tail(As), lattice, hleft, hright)
end

_cu_environments(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, A::GrassmannMPS, B::Vararg{GrassmannMPS}) = cuTwosideExpectationCache(lattice, (A, B...))

