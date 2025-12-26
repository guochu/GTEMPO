
function prodmps(::Type{T}, physpaces::Vector{S}, physectors::Vector; left::S=oneunit(S), right::S=oneunit(S)) where {T <: Number, S <: ElementarySpace}
	L = length(physpaces)
	(L == length(physectors)) || throw(DimensionMismatch())
	physectors = [convert(sectortype(S), item) for item in physectors]

	# the total quantum number is ignored in the Abelian case
	if FusionStyle(sectortype(S)) isa UniqueFusion
		rightind, = ⊗(physectors...)
		right = S((rightind=>1,))
	end
	virtualpaces = Vector{S}(undef, L+1)
	virtualpaces[1] = left
	for i in 2:L
		virtualpaces[i] = fuse(virtualpaces[i-1], S((physectors[i-1]=>1,)) )
	end
	virtualpaces[L+1] = right
	for i in L:-1:2
		virtualpaces[i] = infimum(virtualpaces[i], fuse(virtualpaces[i+1],  S((physectors[i]=>1,))' ))
	end
	return MPS(ones, T, physpaces, virtualpaces)
end
prodmps(::Type{T}, physpace::S, physectors::Vector; kwargs...) where {T <: Number, S <: ElementarySpace} = prodmps(T, [physpace for i in 1:length(physectors)], physectors; kwargs...)
prodmps(physpaces::Vector{S}, physectors::Vector; kwargs...) where {S <: ElementarySpace} = prodmps(Float64, physpaces, physectors; kwargs...)
prodmps(physpace::S, physectors::Vector; kwargs...) where {S <: ElementarySpace} = prodmps(Float64, physpace, physectors; kwargs...)


function prodmpo(::Type{T}, physpaces::Vector{S}, positions::Vector{Int}, ops::Vector) where {T <: Number, S <: ElementarySpace}
	@assert all(x->isa(x, SiteOperator), ops)
	@assert issorted(positions)
	L = length(physpaces)
	for (k, v) in zip(positions, ops)
		((k>= 1) && (k <= L)) || throw(BoundsError(physpaces, k))
		(physpaces[k] == ophysical_space(v) == iphysical_space(v)') || throw(SpaceMismatch("space mismatch on site $k"))
	end
	A = mpotensortype(S, T)
	mpotensors = Vector{A}(undef, L)
	left = oneunit(S)
	for i in 1:L
		pos = findfirst(x->x==i, positions)
		if isnothing(pos)
			mj = id(left ⊗ physpaces[i])
		else
			tmp = ops[pos]
			if isa(tmp, MPSBondTensor)
				mj = _add_legs(tmp, left)
			else
				mj = tmp
			end
		end
		mpotensors[i] = convert(A, mj) 
		left = space_r(mpotensors[i])'
	end
	return MPO(mpotensors)
end





function spin_site_ops_u1()
    ph = Rep[U₁](0=>1, 1=>1)
    vacuum = oneunit(ph)
    σ₊ = zeros(vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
    copy!(block(σ₊, Irrep[U₁](1)), ones(1, 1))
    σ₋ = zeros(vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph)
    copy!(block(σ₋, Irrep[U₁](0)), ones(1, 1))
    σz = ones(ph ← ph)
    copy!(block(σz, Irrep[U₁](0)), -ones(1, 1))
    return Dict("+"=>σ₊, "-"=>σ₋, "z"=>σz)
end
