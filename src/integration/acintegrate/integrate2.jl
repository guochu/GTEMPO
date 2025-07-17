# 2 gmps integration

# function update_pair_left(left::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS; trunc=DefaultIntegrationTruncation)
# 	f = (scaling(x) * scaling(y))^2
# 	posa = 2*j-1
# 	@tensor twositemps1[4,3,5;6] := left[1,2]*x[posa][1,3,4]*y[posa][2,5,6]
# 	mpsj1 = _fuse_physical(swap12!(twositemps1))
# 	m1, m2 = _swap_gate(mpsj1, y[posa+1], trunc=trunc)
# 	@tensor twositemps2[3,2,4;5] := x[posa+1][1,2,3] * m1[1,4,5]
# 	mpsj2 = _fuse_physical(swap12!(twositemps2))
# 	@tensor twositemps3[1,2,4;5] := mpsj2[1,2,3] * m2[3,4,5]
# 	left = _trace_physical(twositemps3, nt=false)
# 	left = rmul!(left, f)
# 	return left	
# end

# function update_pair_right(right::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS; trunc=DefaultIntegrationTruncation)
# 	f = (scaling(x) * scaling(y))^2
# 	posb = 2 * j
# 	@tensor twositemps1[1; 2 6 5] := y[posb][1,2,3] * right[3,4] *x[posb][5,6,4]
# 	mpsj1 = _fuse_physical(swap34!(twositemps1))
# 	m1, m2 = _swap_gate(y[posb-1], mpsj1, trunc=trunc)
# 	@tensor twositemps2[1; 2 5 4] := m2[1,2,3] * x[posb-1][4,5,3] 
# 	mpsj2 = _fuse_physical(swap34!(twositemps2))
# 	@tensor twositemps3[1,2,4;5] := m1[1,2,3] * mpsj2[3,4,5] 
# 	right = _trace_physical(twositemps3, nt=false)
# 	right = rmul!(right, f )
# 	return right
# end

fermionparity(s::Z2Irrep) = isodd(s.n)

function update_pair_left(left::GrassmannTensorMap{<:AbstractTensorMap{<:Number, S, 1, 2}}, j::Int, x::Vector, y::Vector) where {S}
	posa = 2*j-1
	@tensor tmp1[1,4,5;2] := left[1,2,3] * y[posa][3,4,5]
	@tensor tmp2[1,2,5,3;6] := tmp1[1,2,3,4] * x[posa][4,5,6]

	# fuse physical
	tmp3 = g_fuse(tmp2, 2)

	@tensor tmp1[1,2,5,6;3] := tmp3[1,2,3,4] * x[posa+1][4,5,6]
	@tensor tmp2[1,2,3,6;4,7] := tmp1[1,2,3,4,5] * y[posa+1][5,6,7]

	# fuse physical
	tmp3 = g_fuse(tmp2, 3)

	# trace physices
	return permute(g_trace(tmp3, 2), (1,), (2,3))
	# left = TensorMap(zeros, scalartype(tmp3), space(tmp3, 1) ← space(tmp3, 4)' ⊗ space(tmp3, 5)')
	# for (f1, f2) in fusiontrees(tmp3)
	# 	if f1.uncoupled[2] == f1.uncoupled[3]
	# 		f0 = FusionTree((f1.uncoupled[1],), f1.coupled, (f1.isdual[1],))
	# 		@tensor left[f0, f2][1,3,4] += tmp3[f1, f2][1,2,2,3,4]
	# 	end
	# end	

	# return left	
end

function update_pair_right(right::GrassmannTensorMap{<:AbstractTensorMap{<:Number, S, 2, 1}}, j::Int, x::Vector, y::Vector) where S
	posb = 2 * j
	@tensor tmp1[4,1,2;5] := y[posb][1,2,3] * right[3,4,5]
	@tensor tmp2[1,4,2,5; 6] := x[posb][1,2,3] * tmp1[3,4,5,6]

	# fuse physical
	tmp3 = g_fuse(tmp2, 3)

	@tensor tmp1[4,1,2,5;6] := x[posb-1][1,2,3] * tmp3[3,4,5,6]
	@tensor tmp2[1,4,2,5,6;7] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7]

	# fuse physical
	tmp3 = g_fuse(tmp2, 3)

	# trace physices
	return g_trace(tmp3, 3)
	# right = TensorMap(zeros, scalartype(tmp3), space(tmp3, 1) ⊗ space(tmp3, 2) ← space(tmp3, 5)')
	# for (f1, f2) in fusiontrees(tmp3)
	# 	if f1.uncoupled[3] == f1.uncoupled[4]
	# 		f0 = FusionTree((f1.uncoupled[1], f1.uncoupled[2]), f1.coupled, (f1.isdual[1],f1.isdual[2]))
	# 		@tensor right[f0, f2][1,2,4] += tmp3[f1, f2][1,2,3,3,4]
	# 	end
	# end	


	# return right
end

