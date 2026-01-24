"""
	GrassmannMPS{A<:MPSTensor, B<:MPSBondTensor}

Each site of GMPS is divided by a scalar "scaling", such 
that the site tensors of a GMPS is always stored in canonical
form 
It should be understood that each site tensor of the GMPS 
should be multiplied by "scaling"
"""
# struct GrassmannMPS{A <: MPSTensor, B <: MPSBondTensor} <: AbstractFiniteGMPS{A}
# 	data::Vector{A}
# 	svectors::Vector{Union{Missing, B}}
# 	scaling::Ref{Float64}
# end
struct GrassmannMPS{A <: MPSTensor, B <: MPSBondTensor, VA<:AbstractVector{A}} <: AbstractFiniteGMPS{A}
	data::VA
	svectors::Vector{Union{Missing, B}}
	scaling::Ref{Float64}
end

function GrassmannMPS(data::Vector{A}; scaling::Real=1) where {A <: MPSTensor}
	# @assert iseven(length(data))
	isoneunit(space_r(data[end])) || throw(ArgumentError("input data must be in the even sector"))
	T = real(scalartype(A))
	B = bondtensortype(spacetype(A), T)
	svectors = Vector{Union{Missing, B}}(missing, length(data)+1)
	# svectors[1] = Diagonal(id(space_l(data[1])))
	# svectors[end] = Diagonal(id(space_r(data[end])'))
	svectors[1] = DiagonalTensorMap{Float64}(ones, space_l(data[1]) )
	svectors[end] = DiagonalTensorMap{Float64}(ones, space_r(data[end])' )
	return GrassmannMPS(data, svectors, Ref(convert(Float64, scaling)))
end 
function GrassmannMPS(data::AbstractVector{A}, svectors::Vector, scaling::Real=1) where {A <: MPSTensor}
	# @assert iseven(length(data))
	# isoneunit(space_r(data[end])) || throw(ArgumentError("input data must be in the even sector"))
	return GrassmannMPS(data, svectors, Ref(convert(Float64, scaling)))
end 

function GrassmannMPS(::Type{T}, L::Int) where {T <: Number}
	@assert iseven(L)
	s = oneunit(_ph)
	v = zeros(T, s ⊗ _ph ← s )
	copy!(block(v, Z2Irrep(0)), ones(1,1))
	data = [copy(v) for i in 1:L]
	return GrassmannMPS(data, scaling=1)
end
GrassmannMPS(L::Int) = GrassmannMPS(Float64, L)


function Base.getproperty(psi::GrassmannMPS, s::Symbol)
	if s == :s
		return psi.svectors
		# return MPSBondView(psi)
	else
		return getfield(psi, s)
	end
end

# GrassmannMPS(psi::MPS; kwargs...) = GrassmannMPS(psi.data; kwargs...)

"""
	scaling(x::GrassmannMPS)

Return the scaling factor of the GMPS
"""
scaling(x::GrassmannMPS) = x.scaling[]
"""
	setscaling!(x::GrassmannMPS, new_scaling::Real)
Replace the scaling factor of the GMPS by new_scaling 
"""
setscaling!(x::GrassmannMPS, scaling::Real) = (x.scaling[] = scaling)

Base.length(x::GrassmannMPS) = length(x.data)
Base.isempty(x::GrassmannMPS) = isempty(x.data)
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

svectors_uninitialized(psi::GrassmannMPS) = any(ismissing, psi.svectors)
function unset_svectors!(psi::GrassmannMPS)
	psi.svectors[2:end-1] .= missing
	return psi
end

Base.copy(psi::GrassmannMPS) = GrassmannMPS(copy(psi.data), copy(psi.svectors), Ref(scaling(psi)))
function Base.complex(psi::GrassmannMPS)
	if scalartype(psi) <: Real
		data = [complex(item) for item in psi.data]
		return GrassmannMPS(data, psi.svectors, psi.scaling)
	end
	return psi
end

isrightcanonical(a::GrassmannMPS; kwargs...) = (scaling(a) ≈ 1) && all(x->isrightcanonical(x; kwargs...), a.data)
function isleftcanonical(a::GrassmannMPS; kwargs...)
	(scaling(a) ≈ 1) || return false
	all(x->isleftcanonical(x; kwargs...), a.data[1:end-1]) || return false
	return isleftcanonical_r(a.data[end]; kwargs...)
end
function iscanonical(a::GrassmannMPS; kwargs...)
	(scaling(a) ≈ 1) || return false
	isrightcanonical(a) || return false
	# we also check whether the singular vectors are the correct Schmidt numbers
	svectors_uninitialized(a) && return false
	hold = l_LL(a[1], a[1])
	for i in 1:length(a)-1
		hold = updateleft(hold, a[i], a[i])
		tmp = a.s[i+1] * a.s[i+1]
		isapprox(hold, tmp; kwargs...) || return false
	end
	return true
end


# Base.:+(x::GrassmannMPS, y::GrassmannMPS) = GrassmannMPS(x.data + y.data)

# function apply!(t::PartialMPO, mps::GrassmannMPS)
# 	@assert isoneunit(space_r(t))
# 	apply!(t, MPS(mps.data, mps.svectors))
# 	return mps
# 	# _start, _end = positions(t)[1], positions(t)[end]
# 	# S = spacetype(mps)
# 	# T = scalartype(mps)
# 	# physpaces = [space(mps[i], 2) for i in _start:_end]
# 	# mpo = prodmpo(T, physpaces, positions(t) .- (_start-1), op(t)) * TK.scalar(coeff(t))

# 	# M = tensormaptype(S, 2, 3, T)
# 	# r = Vector{M}(undef, _end - _start + 1)
# 	# for (i, pos) in enumerate(_start:_end)
# 	# 	r[i] = @tensor tmp[-1 -2; -3 -4 -5] := mpo[i][-1, -3, -4, 1] * mps[pos][-2, 1, -5]
# 	# end
# 	# fusion_ts = [isomorphism(space(item, 4)' ⊗ space(item, 5)', fuse(space(item, 4)', space(item, 5)')) for item in r]
# 	# left = isomorphism(fuse(space(r[1], 1), space(r[1], 2)), space(r[1], 1) ⊗ space(r[1], 2))
# 	# mps[_start] = @tensor tmp[1,4;7] := left[1,2,3] * r[1][2,3,4,5,6] * fusion_ts[1][5,6,7]
# 	# for (i, pos) in enumerate(_start+1:_end)
# 	# 	mps[pos] = @tensor tmp[3,4;7] := conj(fusion_ts[i][1,2,3]) * r[i+1][1,2,4,5,6] * fusion_ts[i+1][5,6,7]
# 	# end
# 	# return mps
# end
function apply!(t::PartialMPO, mps::GrassmannMPS)
	@assert isoneunit(space_r(t))
    @assert positions(t)[end] <= length(mps)
    T = promote_type(scalartype(t), scalartype(mps))
    M = tensormaptype(spacetype(mps), 2, 3, T)
    _start, _end = positions(t)[1], positions(t)[end]
    r = Vector{M}(undef, _end - _start + 1)
    leftspace = oneunit(space_l(t))

    for (i, pos) in enumerate(_start:_end)
        _pos = findfirst(x->x==pos, positions(t))
        if !isnothing(_pos)
            r[i] = @tensor tmp[-1 -2; -3 -4 -5] := t[_pos][-1, -3, -4, 1] * mps[pos][-2, 1, -5]
            leftspace = space_r(t[_pos])'
        else
            hj = id(storagetype(M), leftspace)
            r[i] = @tensor tmp[-1 -2; -3 -4 -5] := hj[-1, -4] * mps[pos][-2, -3, -5]
        end
    end
    fusion_ts = [isomorphism(T, space(item, 4)' ⊗ space(item, 5)', fuse(space(item, 4)', space(item, 5)')) for item in r]
    left = isomorphism(T, fuse(space(r[1], 1), space(r[1], 2)), space(r[1], 1) ⊗ space(r[1], 2))
    mps[_start] = @tensor tmp[1,4;7] := left[1,2,3] * r[1][2,3,4,5,6] * fusion_ts[1][5,6,7]
    for (i, pos) in enumerate(_start+1:_end)
        mps[pos] = @tensor tmp[3,4;7] := conj(fusion_ts[i][1,2,3]) * r[i+1][1,2,4,5,6] * fusion_ts[i+1][5,6,7]
    end
    unset_svectors!(mps)
    return mps
end
Base.:*(t::PartialMPO, mps::GrassmannMPS) = apply!(t, copy(mps))

function apply!(x::Union{GTerm, ExpGTerm}, mps::GrassmannMPS)
	all(s -> 1 <= s <= length(mps), positions(x)) || throw(BoundsError())
	t = convert(PartialMPO, x)
	apply!(t, mps)
	return mps
end
Base.:*(x::Union{GTerm, ExpGTerm}, mps::GrassmannMPS) = apply!(x, copy(mps))
function apply!(x::GTerm{0}, mps::GrassmannMPS)
	return mps * x.coeff
end

# function randomgmps(::Type{T}, L::Int; D::Int) where {T <: Number}
# 	physpaces = [_ph for i in 1:L]
# 	mps = randommps(T, physpaces, D=D)
# 	return GrassmannMPS(mps.data)
# end

function randomgmps(::Type{T}, L::Int; D::Int) where {T <: Number}
	Dh = div(D, 2)
	virtualspace = Z2Space(0=>Dh, 1=>Dh)
	virtualspaces = Vector{typeof(virtualspace)}(undef, L+1)
	virtualspaces[1] = virtualspaces[end] = Z2Space(0=>1)
	for i in 2:L
		virtualspaces[i] = virtualspace
	end
	physpaces = [_ph for i in 1:L]

	# data = [TensorMap(randn, T, virtualspaces[i] ⊗ physpaces[i] ← virtualspaces[i+1]) for i in 1:L]
	data = [randn(T, virtualspaces[i] ⊗ physpaces[i] ← virtualspaces[i+1]) for i in 1:L]
	for i in 1:L
		data[i] /= sqrt(D)
	end
	return GrassmannMPS(data)
end
randomgmps(L::Int; kwargs...) = randomgmps(Float64, L; kwargs...)


function increase_bond!(x::GrassmannMPS, D::Int)
	Dh = div(D, 2)
	virtualspace = Z2Space(0=>Dh, 1=>Dh)
	L = length(x)
	ms = [isometry(scalartype(x), virtualspace, space_l(x[site+1])) for site in 1:L-1]
	@tensor tmp[1,2;4] := x[1][1,2,3] * conj(ms[1][4,3])
	x[1] = tmp
	@tensor tmp[1,3;4] := ms[L-1][1,2] * x[L][2,3,4]
	x[L] = tmp
	for site in 2:L-1
		vspace = space_l(x[site+1])
		(vspace ≾ virtualspace) || @warn "old virtualspace $(vspace) is not monomorphic to the new virtualspace $(virtualspace)"
		@tensor tmp[1,3;5] := ms[site-1][1,2] * x[site][2,3,4] * conj(ms[site][5,4])
		x[site] = tmp
	end
	unset_svectors!(x)
	return x
end

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
# 	u, s, v, err = stable_tsvd!(twositemps1; trunc=trunc)
# 	x.data.s[bond+1] = s
# 	x[bond] = twositemps * v' 
# 	x[bond+1] = permute(v, (1,2), (3,))
# 	return x
# end

function easy_swap!(x::AbstractGMPS, bond::Int; trunc::TruncationScheme=DefaultTruncation)
	x[bond], x.s[bond+1], x[bond+1] = _swap_gate(x.s[bond], x[bond], x.s[bond+1], x[bond+1], trunc=trunc)
	# x[bond], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
	return x
end

function naive_swap!(x::AbstractGMPS, bond::Int; trunc::TruncationScheme=DefaultTruncation)
	x[bond], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
	# x[bond], x[bond+1] = _swap_gate(x[bond], x[bond+1], trunc=trunc)
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

# function _fuse_physical(m::AbstractTensorMap{S, 4, 2}) where S
# 	cod = ProductSpace{S}((space(m, 1), space(m, 2), space(m, 3)))
# 	dom = ProductSpace{S}((space(m, 5)', space(m, 6)'))
# 	r = TensorMap(ds->zeros(scalartype(m), ds), cod ← dom)
# 	for (f1, f2) in fusiontrees(m)
# 		if f1.uncoupled[3] == f1.uncoupled[4]
# 			if iseven(f1.uncoupled[3].n)
# 				f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[3]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3]))
# 				r[f1′, f2] .+= m[f1, f2][:, :, :, 1, :, :]
# 			end
# 		else
# 			isdual3 = isodd(f1.uncoupled[3].n) ? f1.isdual[3] : f1.isdual[4]
# 			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],Z2Irrep(1)), f1.coupled, (f1.isdual[1],f1.isdual[2], isdual3))
# 			r[f1′, f2] .+= m[f1, f2][:, :, :, 1, :, :]
# 		end
# 	end
# 	return r
# end

function _swap_gate(m1, m2; trunc)
	@tensor twositemps[1,4;2,5] := GrassmannTensorMap(m1)[1,2,3] * GrassmannTensorMap(m2)[3,4,5]
	u, s, v, err = stable_tsvd!(twositemps; trunc=trunc)
	# return u, permute(s * v, (1,2), (3,))
	return get_data(u * s), get_data(permute(v, (1,2), (3,)))
end

function _swap_gate(svectorj1, m1, svectorj2, m2; trunc::TruncationScheme)
	@tensor twositemps[1,4;2,5] := GrassmannTensorMap(m1)[1,2,3] * GrassmannTensorMap(m2)[3,4,5]
	# println(space(svectorj1, 2), " ", space(m1, 1))
	@tensor twositemps1[-1 -2; -3 -4] := GrassmannTensorMap(svectorj1)[-1, 1] * twositemps[1, -2, -3, -4]
	u, s, v, err = stable_tsvd!(twositemps1, trunc=trunc)
	@tensor u[-1 -2; -3] = twositemps[-1,-2,1,2] * conj(v[-3,1,2])
	return get_data(u), get_data(s), get_data(permute(v, (1,2), (3,)))
end

