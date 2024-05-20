# integrate multiple GMPSs using BMPS-like method

function _ac_integrate(alg::BMPSIntegrate, lattice::AbstractGrassmannLattice, x0::GrassmannMPS, x1::GrassmannMPS, x2::GrassmannMPS...; center::Int=pos2pairindex(length(lattice)))
	(ConjugationStyle(lattice) isa AdjacentConjugation) || throw(ArgumentError("AdjacentConjugation expected"))
	xx = collect((x0, x1, x2...))
	return _bmps_integrate(lattice, xx; center=center, trunc=alg.trunc)
end

function _bmps_integrate(lattice::AbstractGrassmannLattice, xx::Vector{<:GrassmannMPS}; center::Int=pos2pairindex(length(lattice)), kwargs...)
	left = _bmps_integrate_left( xx, center+1; kwargs...)
	right = _bmps_integrate_right(xx, center; kwargs...)
	return contract_center(left, right)
end

function _integrate_center(xx::Vector{<:GrassmannMPS}, left::GrassmannMPS, left_b::Int, 
		right::GrassmannMPS, right_b::Int, center::Int = div(left_b+right_b, 2); trunc::TruncationScheme=DefaultIntegrationTruncation)
	@assert left_b <= center <= right_b
	for i in left_b:center
		left = update_pair_left(left, i, xx; trunc=trunc)
	end
	for i in right_b:-1:center+1
		right = update_pair_right(right, i, xx; trunc=trunc)
	end
	return contract_center(left, right)
end

function _bmps_integrate_left(xx::Vector{<:GrassmannMPS}, center::Int; trunc::TruncationScheme=DefaultIntegrationTruncation)
	left = _l_bmps_boundary(xx)
	for i in 1:center-1
		left = update_pair_left(left, i, xx, trunc=trunc)
	end
	return left
end
function _bmps_integrate_right(xx::Vector{<:GrassmannMPS}, center::Int; trunc::TruncationScheme=DefaultIntegrationTruncation)
	right = _r_bmps_boundary(xx)
	Lhalf = pos2pairindex(length(xx[1]))
	for i in Lhalf:-1:center+1
		right = update_pair_right(right, i, xx, trunc=trunc)
	end	
	return right
end

function _l_bmps_boundary(xx::Vector{<:GrassmannMPS})
	s = oneunit(_ph)
	v = TensorMap(ds->ones(scalartype(xx[1]), ds), s ⊗ s' ← s )
	return GrassmannMPS([v for i in 1:length(xx)])
end

function _r_bmps_boundary(xx::Vector{<:GrassmannMPS})
	s = oneunit(_ph)
	v = TensorMap(ds->ones(scalartype(xx[1]), ds), s ⊗ s ← s )
	return GrassmannMPS([v for i in 1:length(xx)])
end

function update_left_util(left::GrassmannMPS, pos::Int, x::Vector{<:GrassmannMPS}; trunc=DefaultIntegrationTruncation)
	@assert length(left) == length(x)
	tmp = [rmul!(_apply_physical_left(left[i], x[i][pos]), scaling(x[i]) ) for i in 1:length(left)]
	left2 = similar(left.data, length(left))
	for i in length(left)-1:-1:1
		a, b = swap_left(tmp[i], tmp[i+1], trunc)
		left2[i+1] = get_data(b)
		tmp[i] = a
	end
	return tmp, left2
end

function update_left_1(left::GrassmannMPS, pos::Int, x::Vector{<:GrassmannMPS}; trunc)
	tmp, left2 = update_left_util(left, pos, x)
	# fuse boundary
	left2[1] = _fuse_boundary(tmp[1])
	return canonicalize!(GrassmannMPS(left2, scaling=scaling(left)), alg=Orthogonalize(trunc=trunc))
end
function update_left_2(left::GrassmannMPS, pos::Int, x::Vector{<:GrassmannMPS}; trunc)
	tmp, left2 = update_left_util(left, pos, x)
	# trace physices
	left2[1] = _trace_boundary(tmp[1])
	return canonicalize!(GrassmannMPS(left2, scaling=scaling(left)), alg=Orthogonalize(trunc=trunc))
end

function update_pair_left(left::GrassmannMPS, j::Int, xx::Vector{<:GrassmannMPS}; trunc=DefaultIntegrationTruncation)
	posa = 2 * j -1
	left = update_left_1(left, posa, xx, trunc=trunc)
	left = update_left_2(left, posa+1, xx, trunc=trunc)
	return left
end

function update_right_util(right::GrassmannMPS, pos::Int, x::Vector{<:GrassmannMPS}; trunc=DefaultIntegrationTruncation)
	@assert length(right) == length(x)
	L = length(right)
	tmp = [rmul!(_apply_physical_right(right[i], x[L-i+1][pos]), scaling(x[L-i+1])) for i in 1:L]
	right2 = similar(right.data, L)
	for i in length(right)-1:-1:1
		a, b = swap_left(tmp[i], tmp[i+1], trunc)
		right2[i+1] = get_data(b)
		tmp[i] = a
	end
	return tmp, right2
end

function update_right_1(right::GrassmannMPS, pos::Int, x::Vector{<:GrassmannMPS}; trunc)
	tmp, right2 = update_right_util(right, pos, x)
	right2[1] = _fuse_boundary(tmp[1])
	return canonicalize!(GrassmannMPS(right2, scaling=scaling(right)), alg=Orthogonalize(trunc=trunc))
end
function update_right_2(right::GrassmannMPS, pos::Int, x::Vector{<:GrassmannMPS}; trunc)
	tmp, right2 = update_right_util(right, pos, x)
	right2[1] = _trace_boundary(tmp[1], nt=false)
	return canonicalize!(GrassmannMPS(right2, scaling=scaling(right)), alg=Orthogonalize(trunc=trunc))
end

function update_pair_right(right::GrassmannMPS, j::Int, xx::Vector{<:GrassmannMPS}; trunc)
	posb = 2 * j
	right = update_right_1(right, posb, xx, trunc=trunc)
	right = update_right_2(right, posb-1, xx, trunc=trunc)
	return right	
end

function contract_center(left::GrassmannMPS, right::GrassmannMPS)
	L = length(left)
	mj = GrassmannTensorMap(isomorphism(Matrix{scalartype(left)}, space_r(left)', space_l(right)))
	f = scaling(left) * scaling(right)
	for i in L:-1:1
		mj = @tensor tmp[1,5] := f * GrassmannTensorMap(left[i])[1,2,3] * mj[3,4] * GrassmannTensorMap(right[L-i+1])[4,2,5]
	end
	return TK.scalar(mj)
end

function _fuse_boundary(tmp::GrassmannTensorMap{<:AbstractTensorMap{S, 3, 1}}) where S
	tmp1 = tmp.data
	m = isomorphism(Matrix{scalartype(tmp1)}, fuse(space(tmp1, 1), space(tmp1, 2)), space(tmp1, 1) ⊗ space(tmp1, 2))
	@tensor l[1,4;5] := m[1,2,3] * tmp1[2,3,4,5]
	return l	
end

function _trace_boundary(tmp1::GrassmannTensorMap{<:AbstractTensorMap{S, 3, 1}}; nt::Bool=true) where S
	# trace physices
	tmp1 = get_data(tmp1)
	m1 = TensorMap(ds->zeros(scalartype(tmp1), ds), space(tmp1, 3) ← space(tmp1, 4)' )
	for (f1, f2) in fusiontrees(tmp1)
		if f1.uncoupled[1] == f1.uncoupled[2]
			f0 = FusionTree((f1.uncoupled[3],), f1.coupled, (f1.isdual[3],))
			coef = (isodd(f1.uncoupled[1].n) && (!nt)) ? -1 : 1
			@tensor m1[f0, f2][2,3] += coef * tmp1[f1, f2][1,1,2,3]
		end
	end	
	vacuum = oneunit(spacetype(tmp1))
	util = TensorMap(ones, scalartype(tmp1), vacuum)
	@tensor l[1,2;3] := util[1] * m1[2,3]
	return l
end

function _apply_physical_left(a::MPSTensor, b::MPSTensor)
	# println(space(a1, 2), " ", space(b, 1))
	@tensor c[1,4,5;3] := GrassmannTensorMap(a)[1,2,3] * GrassmannTensorMap(b)[2,4,5]
	return c
end

function _apply_physical_right(a::MPSTensor, b::MPSTensor)
	@tensor c[4,2,1;5] := GrassmannTensorMap(b)[1,2,3] * GrassmannTensorMap(a)[4,3,5]
	return c
end

function swap_left(a::GrassmannTensorMap{<:AbstractTensorMap{S, 3, 1}}, b::GrassmannTensorMap{<:AbstractTensorMap{S, 3, 1}}, trunc) where S
	@tensor tmp2[1,2,5,3;6,7] := a[1,2,3,4] * b[4,5,6,7]
	# fuse indices
	# cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) 
	# dom = space(tmp2, 5)' ⊗ space(tmp2, 6)'
	# tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	# for (f1, f2) in fusiontrees(tmp2)
	# 	n = f1.uncoupled[2].n + f1.uncoupled[3].n
	# 	if n == 0
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[4]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[4]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:]
	# 	elseif n == 1
	# 		isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
	# 		f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4]), f1.coupled, (f1.isdual[1], isdual2, f1.isdual[4]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:]
	# 	end
	# end
	tmp3 = g_fuse(tmp2, 2)
	u, s, v, err = stable_tsvd!(tmp3; trunc=trunc)
	return u, permute(s * v, (1,2), (3,))
end
