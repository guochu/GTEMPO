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

function update_pair_left(left::AbstractTensorMap{<:ElementarySpace, 1, 2}, j::Int, x::Vector, y::Vector)
	posa = 2*j-1
	@tensor tmp1[1,4,5;2] := left[1,2,3] * y[posa][3,4,5]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (fermionparity(f1.uncoupled[2]) && fermionparity(f2.uncoupled[1])) ? -1 : 1
		coef2 = (fermionparity(f1.uncoupled[3]) && fermionparity(f2.uncoupled[1])) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	
	@tensor tmp2[1,2,5,3;6] := tmp1[1,2,3,4] * x[posa][4,5,6]
	for (f1, f2) in fusiontrees(tmp2)
		coef = (fermionparity(f1.uncoupled[3]) && fermionparity(f1.uncoupled[4])) ? -1 : 1
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	# fuse physical
	# cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4)
	# dom = space(tmp2, 5)'
	# tmp3 = TensorMap(zeros, scalartype(tmp2), cod ← dom) 
	# for (f1, f2) in fusiontrees(tmp2)
	# 	n = f1.uncoupled[2].n + f1.uncoupled[3].n 
	# 	if n == 0
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[4]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[4]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:]
	# 	elseif n == 1
	# 		isdual2 = ifelse(fermionparity(f1.uncoupled[2]), f1.isdual[2], f1.isdual[3] )
	# 		f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4]), f1.coupled, (f1.isdual[1], isdual2, f1.isdual[4]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:]
	# 	end
	# end
	tmp3 = g_fuse(tmp2, 2)

	@tensor tmp1[1,2,5,6;3] := tmp3[1,2,3,4] * x[posa+1][4,5,6]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (fermionparity(f1.uncoupled[3]) && fermionparity(f2.uncoupled[1])) ? -1 : 1
		coef2 = (fermionparity(f1.uncoupled[4]) && fermionparity(f2.uncoupled[1])) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	

	@tensor tmp2[1,2,3,6;4,7] := tmp1[1,2,3,4,5] * y[posa+1][5,6,7]
	for (f1, f2) in fusiontrees(tmp2)
		coef = (fermionparity(f1.uncoupled[4]) && fermionparity(f2.uncoupled[1])) ? -1 : 1
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	# fuse physical
	# cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) 
	# dom = space(tmp2, 5)' ⊗ space(tmp2, 6)'
	# tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	# for (f1, f2) in fusiontrees(tmp2)
	# 	n = f1.uncoupled[3].n + f1.uncoupled[4].n 
	# 	if n == 0
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:]
	# 	elseif n == 1
	# 		isdual2 = ifelse(fermionparity(f1.uncoupled[3]), f1.isdual[3], f1.isdual[4] )
	# 		f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2], Z2Irrep(1)), f1.coupled, (f1.isdual[1], f1.isdual[2], isdual2))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:]
	# 	end
	# end
	tmp3 = g_fuse(tmp2, 3)

	# trace physices
	left = TensorMap(zeros, scalartype(tmp3), space(tmp3, 1) ← space(tmp3, 4)' ⊗ space(tmp3, 5)')
	for (f1, f2) in fusiontrees(tmp3)
		if f1.uncoupled[2] == f1.uncoupled[3]
			f0 = FusionTree((f1.uncoupled[1],), f1.coupled, (f1.isdual[1],))
			@tensor left[f0, f2][1,3,4] += tmp3[f1, f2][1,2,2,3,4]
		end
	end	

	return left	
end

function update_pair_right(right::AbstractTensorMap{<:ElementarySpace, 2, 1}, j::Int, x::Vector, y::Vector)
	posb = 2 * j
	@tensor tmp1[4,1,2;5] := y[posb][1,2,3] * right[3,4,5]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (fermionparity(f1.uncoupled[1]) && fermionparity(f1.uncoupled[2])) ? -1 : 1
		coef2 = (fermionparity(f1.uncoupled[1]) && fermionparity(f1.uncoupled[3])) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	
	@tensor tmp2[1,4,2,5; 6] := x[posb][1,2,3] * tmp1[3,4,5,6]
	for (f1, f2) in fusiontrees(tmp2)
		coef = (fermionparity(f1.uncoupled[2]) && fermionparity(f1.uncoupled[3])) ? -1 : 1
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	# fuse physical
	# cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) 
	# dom = space(tmp2, 5)' 
	# tmp3 = TensorMap(zeros, scalartype(tmp2), cod ← dom) 
	# for (f1, f2) in fusiontrees(tmp2)
	# 	n = f1.uncoupled[3].n + f1.uncoupled[4].n 
	# 	if n == 0
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:]
	# 	elseif n == 1
	# 		isdual2 = ifelse(fermionparity(f1.uncoupled[3]), f1.isdual[3], f1.isdual[4] )
	# 		f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2], Z2Irrep(1)), f1.coupled, (f1.isdual[1], f1.isdual[2], isdual2))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:]
	# 	end
	# end
	tmp3 = g_fuse(tmp2, 3)

	@tensor tmp1[4,1,2,5;6] := x[posb-1][1,2,3] * tmp3[3,4,5,6]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (fermionparity(f1.uncoupled[1]) && fermionparity(f1.uncoupled[2])) ? -1 : 1
		coef2 = (fermionparity(f1.uncoupled[1]) && fermionparity(f1.uncoupled[3])) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	

	@tensor tmp2[1,4,2,5,6;7] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7]
	for (f1, f2) in fusiontrees(tmp2)
		coef = (fermionparity(f1.uncoupled[2]) && fermionparity(f1.uncoupled[3])) ? -1 : 1
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	# fuse physical
	# cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 5) 
	# dom = space(tmp2, 6)' 
	# tmp3 = TensorMap(zeros, scalartype(tmp2), cod ← dom) 
	# for (f1, f2) in fusiontrees(tmp2)
	# 	n = f1.uncoupled[3].n + f1.uncoupled[4].n 
	# 	if n == 0
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[5]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3], f1.isdual[5]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:]
	# 	elseif n == 1
	# 		isdual2 = ifelse(fermionparity(f1.uncoupled[3]), f1.isdual[3], f1.isdual[4] )
	# 		f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5]), f1.coupled, (f1.isdual[1], f1.isdual[2], isdual2, f1.isdual[5]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:]
	# 	end
	# end
	tmp3 = g_fuse(tmp2, 3)

	# trace physices
	right = TensorMap(zeros, scalartype(tmp3), space(tmp3, 1) ⊗ space(tmp3, 2) ← space(tmp3, 5)')
	for (f1, f2) in fusiontrees(tmp3)
		if f1.uncoupled[3] == f1.uncoupled[4]
			f0 = FusionTree((f1.uncoupled[1], f1.uncoupled[2]), f1.coupled, (f1.isdual[1],f1.isdual[2]))
			@tensor right[f0, f2][1,2,4] += tmp3[f1, f2][1,2,3,3,4]
		end
	end	


	return right
end

