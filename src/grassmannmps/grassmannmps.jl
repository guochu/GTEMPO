"""
	GrassmannMPS{A, B}
"""
struct GrassmannMPS{A <: MPSTensor} 
	data::Vector{A}
	scale::Ref{Float64}
end

function GrassmannMPS(::Type{T}, L::Int) where {T <: Number}
	@assert iseven(L)
	s = oneunit(_ph)
	v = TensorMap(ds->zeros(T, ds), s ⊗ _ph ← s )
	blocks(v)[Irrep[ℤ₂](0)] = ones(1,1)
	return GrassmannMPS([copy(v) for i in 1:L], Ref(1.))
end
GrassmannMPS(L::Int) = GrassmannMPS(Float64, L)
GrassmannMPS(data::Vector{<:MPSTensor}; scale::Real=1) = GrassmannMPS(data, Ref(convert(Float64, scale)))
# GrassmannMPS(psi::MPS; kwargs...) = GrassmannMPS(psi.data; kwargs...)
scale(x::GrassmannMPS) = x.scale[]

Base.length(x::GrassmannMPS) = length(x.data)
Base.isempty(x::GrassmannMPS) = isempty(x.data)
TK.scalartype(::Type{GrassmannMPS{A}}) where {A} = scalartype(A)
Base.getindex(x::GrassmannMPS, i::Int) = getindex(x.data, i)
Base.firstindex(x::GrassmannMPS) = firstindex(x.data)
Base.lastindex(x::GrassmannMPS) = lastindex(x.data)
function Base.setindex!(x::GrassmannMPS, v::MPSTensor, i::Int)
	# DMRG.check_mpstensor_dir(v) || throw(SpaceMismatch())
	# if i == 1
	# 	isoneunit(space_l(v)) || throw(SpaceMismatch("space_l of MPS should be vacuum by convention."))
	# end
	return setindex!(x.data, v, i)
end
Base.copy(x::GrassmannMPS) = GrassmannMPS(copy(x.data), scale=scale(x))
DMRG.bond_dimension(a::GrassmannMPS, bond::Int) = begin
	((bond >= 1) && (bond <= length(a))) || throw(BoundsError())
	dim(space(a[bond], 3))
end 
DMRG.bond_dimensions(a::GrassmannMPS) = [bond_dimension(a, i) for i in 1:length(a)]
DMRG.bond_dimension(a::GrassmannMPS) = maximum(bond_dimensions(a))
TK.spacetype(::Type{GrassmannMPS{A}}) where A = spacetype(A)
TK.spacetype(x::GrassmannMPS) = spacetype(typeof(x))

# Base.:+(x::GrassmannMPS, y::GrassmannMPS) = GrassmannMPS(x.data + y.data)

function apply!(t::QTerm, mps::GrassmannMPS)
	_start, _end = positions(t)[1], positions(t)[end]
	S = spacetype(mps)
	T = scalartype(mps)
	physpaces = [space(mps[i], 2) for i in _start:_end]
	mpo = prodmpo(T, physpaces, positions(t) .- (_start-1), op(t)) * TK.scalar(coeff(t))

	M = tensormaptype(S, 2, 3, T)
	r = Vector{M}(undef, _end - _start + 1)
	for (i, pos) in enumerate(_start:_end)
		r[i] = @tensor tmp[-1 -2; -3 -4 -5] := mpo[i][-1, -3, -4, 1] * mps[pos][-2, 1, -5]
	end
	fusion_ts = [isomorphism(space(item, 4)' ⊗ space(item, 5)', fuse(space(item, 4)', space(item, 5)')) for item in r]
	left = isomorphism(fuse(space(r[1], 1), space(r[1], 2)), space(r[1], 1) ⊗ space(r[1], 2))
	mps[_start] = @tensor tmp[1,4;7] := left[1,2,3] * r[1][2,3,4,5,6] * fusion_ts[1][5,6,7]
	for (i, pos) in enumerate(_start+1:_end)
		mps[pos] = @tensor tmp[3,4;7] := conj(fusion_ts[i][1,2,3]) * r[i+1][1,2,4,5,6] * fusion_ts[i+1][5,6,7]
	end
	return mps
end
Base.:*(t::QTerm, mps::GrassmannMPS) = apply!(t, copy(mps))

function apply!(x::Union{GTerm, ExpGTerm}, mps::GrassmannMPS)
	all(s -> 1 <= s <= length(mps), positions(x)) || throw(BoundsError())
	t = convert(QTerm, x)
	apply!(t, mps)
	return mps
end
Base.:*(x::Union{GTerm, ExpGTerm}, mps::GrassmannMPS) = apply!(x, copy(mps))

# function randomgmps(::Type{T}, L::Int; D::Int) where {T <: Number}
# 	physpaces = [_ph for i in 1:L]
# 	mps = randommps(T, physpaces, D=D)
# 	return GrassmannMPS(mps.data)
# end

function randomgmps(::Type{T}, L::Int; D::Int) where {T <: Number}
	Dh = div(D, 2)
	virtualspace = Rep[ℤ₂](0=>Dh, 1=>Dh)
	virtualspaces = Vector{typeof(virtualspace)}(undef, L+1)
	virtualspaces[1] = virtualspaces[end] = Rep[ℤ₂](0=>1)
	for i in 2:L
		virtualspaces[i] = virtualspace
	end
	physpaces = [_ph for i in 1:L]
	mps = MPS(randn, T, physpaces, virtualspaces)
	for i in 1:L
		mps[i] /= sqrt(D)
	end
	# canonicalize!(mps, alg=Orthogonalize(QR(), normalize=true))
	return GrassmannMPS(mps.data)
	# out = GrassmannMPS(mps.data)
	# return canonicalize!(out, alg=Orthogonalize(QR()))
end
randomgmps(L::Int; kwargs...) = randomgmps(Float64, L; kwargs...)

# # swap the i-th bond
# function swap!(x::GrassmannMPS, bond::Int; trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	# @assert space(x.data.s[bond], 2)' == space_l(x[bond]) 
# 	@assert (1 <= bond < length(x))
# 	@tensor twositemps[1,4;2,5] := x[bond][1,2,3] * x[bond+1][3,4,5]
# 	for (f1, f2) in fusiontrees(twositemps)
# 		if isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)
# 			twositemps[f1, f2] .*= -1
# 		end
# 	end
# 	@tensor twositemps1[1,3;4,5] := x.data.s[bond][1,2] * twositemps[2,3,4,5]
# 	u, s, v, err = DMRG.stable_tsvd!(twositemps1; trunc=trunc)
# 	x.data.s[bond+1] = s
# 	x[bond] = twositemps * v' 
# 	x[bond+1] = permute(v, (1,2), (3,))
# 	return x
# end

function easy_swap!(x::GrassmannMPS, bond::Int; trunc::TruncationScheme=DMRG.DefaultTruncation)
	x[bond], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
	return x
end

# function _mult_site(xj, yj)
# 	S = spacetype(xj)
# 	zj = permute(xj, (1,2,3), ()) ⊗ permute(yj, (1,2,3), ())
# 	p1 = (1,4,2,5)
# 	p2 = (3,6)
#     cod = ProductSpace{S}(map(n->space(zj, n), p1))
#     dom = ProductSpace{S}(map(n->dual(space(zj, n)), p2))
# 	return TK._add!(true, zj, false, similar(zj, cod←dom), p1, p2, (f1, f2)->f_permute(f1, f2, p1, p2))
# end
# # f2 is empty
# function f_permute(f1, f2, p1, p2)
# 	ft, coeff = first(permute(f1, f2, p1, p2))
# 	t4 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
# 	return TK.SingletonDict(ft=>coeff*t4)
# end
function _mult_site(xj, yj)
	@tensor r[1,4,2,5;3,6] := xj[1,2,3] * yj[4,5,6]
	for (f1, f2) in fusiontrees(r)
		if isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)
			r[f1, f2] .*= -1
		end
	end
	return r
end
function _fuse_physical(m::AbstractTensorMap{S, 4, 2}) where S
	cod = ProductSpace{S}((space(m, 1), space(m, 2), space(m, 3)))
	dom = ProductSpace{S}((space(m, 5)', space(m, 6)'))
	r = TensorMap(ds->zeros(scalartype(m), ds), cod ← dom)
	for (f1, f2) in fusiontrees(m)
		if f1.uncoupled[3] == f1.uncoupled[4]
			if iseven(f1.uncoupled[3].n)
				f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[3]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3]))
				r[f1′, f2] .+= m[f1, f2][:, :, :, 1, :, :]
			end
		else
			isdual3 = isodd(f1.uncoupled[3].n) ? f1.isdual[3] : f1.isdual[4]
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],Z2Irrep(1)), f1.coupled, (f1.isdual[1],f1.isdual[2], isdual3))
			r[f1′, f2] .+= m[f1, f2][:, :, :, 1, :, :]
		end
	end
	return r
end

function _swap_gate(m1, m2; trunc)
	@tensor twositemps[1,4;2,5] := m1[1,2,3] * m2[3,4,5]
	for (f1, f2) in fusiontrees(twositemps)
		if isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)
			twositemps[f1, f2] .*= -1
		end
	end
	u, s, v, err = DMRG.stable_tsvd!(twositemps; trunc=trunc)
	# return u, permute(s * v, (1,2), (3,))
	return u * s, permute(v, (1,2), (3,))
end
