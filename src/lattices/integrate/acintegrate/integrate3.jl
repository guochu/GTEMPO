
function update_pair_left(left0::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS; trunc)
	posa = 2*j-1
	f = (scaling(x) * scaling(y) * scaling(z) )^2
	left = copy(left0)

	for (f1, f2) in fusiontrees(left)
		coef1 = (isodd(f1.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef = coef1 * coef2
		if coef != 1
			left[f1, f2] .*= coef
		end
	end
	@tensor tmp1[4,5,3;2] := left[1,2,3] * x[posa][1,4,5] 
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		# println(coef1, " ", coef2, " ", coef3, " ", coef4, " ", coef5)
		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end
	@tensor tmp2[1,5,2,6;3] := tmp1[1,2,3,4] * y[posa][4,5,6]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end		
	@tensor threesitemps1[1,2,6,3,4;7] := tmp2[1,2,3,4,5] * z[posa][5,6,7]
	for (f1, f2) in fusiontrees(threesitemps1)
		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			threesitemps1[f1, f2] .*= coef
		end			
	end

	# fuse indices
	cod = space(threesitemps1, 1) ⊗ space(threesitemps1, 4) ⊗ space(threesitemps1, 5)
	dom = space(threesitemps1, 6)'
	left4 = TensorMap(ds->zeros(scalartype(threesitemps1), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(threesitemps1)
		n = f1.uncoupled[1].n + f1.uncoupled[2].n + f1.uncoupled[3].n
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, (f1.isdual[1],f1.isdual[4],f1.isdual[5]))
			left4[f1′, f2] .+= threesitemps1[f1, f2][:,1,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[1].n), f1.isdual[1], ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3] ))
			f1′ = FusionTree((Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, (isdual2, f1.isdual[4], f1.isdual[5]))
			left4[f1′, f2] .+= threesitemps1[f1, f2][:,1,1,:,:,:]
		end
	end

	# \bar{a}

	for (f1, f2) in fusiontrees(left4)
		coef1 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			left4[f1, f2] .*= coef
		end
	end
	@tensor tmp1[1,5,6,4;3] := left4[1,2,3,4] * x[posa+1][2,5,6]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	
	@tensor tmp2[1,2,6,3,7;4] := tmp1[1,2,3,4,5] * y[posa+1][5,6,7]	
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	@tensor threesitemps2[1,2,3,7,4,5;8] := tmp2[1,2,3,4,5,6] * z[posa+1][6,7,8]
	for (f1, f2) in fusiontrees(threesitemps2)
		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			threesitemps2[f1, f2] .*= coef
		end
	end

	# fuse indices
	cod = space(threesitemps2, 1) ⊗ space(threesitemps2, 2) ⊗ space(threesitemps2, 5) ⊗ space(threesitemps2, 6)
	dom = space(threesitemps2, 7)'
	left5 = TensorMap(ds->zeros(scalartype(threesitemps2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(threesitemps2)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n + f1.uncoupled[4].n
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2],f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[5],f1.isdual[6]))
			left5[f1′, f2] .+= threesitemps2[f1, f2][:,:,1,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4] ))
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, (f1.isdual[1], isdual2, f1.isdual[5], f1.isdual[6]))
			left5[f1′, f2] .+= threesitemps2[f1, f2][:,:,1,1,:,:,:]
		end
	end
	# trace physices
	left = TensorMap(ds->zeros(scalartype(left5), ds), space(left5, 3) ⊗ space(left5, 4) ← space(left5, 5)' )
	for (f1, f2) in fusiontrees(left5)
		if f1.uncoupled[1] == f1.uncoupled[2]
			f0 = FusionTree((f1.uncoupled[3],f1.uncoupled[4]), f1.coupled, (f1.isdual[3],f1.isdual[4]))
			@tensor left[f0, f2][2,3,4] += left5[f1, f2][1,1,2,3,4]
			# axpy!(1, tmp, r[f0, f2]) 
		end
	end	
	left = rmul!(left, f)	

	return left	
end

function update_pair_right(right::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS; trunc)
	posb = 2 * j
	f = (scaling(x) * scaling(y) * scaling(z))^2

	@tensor tmp1[4 ;1 2 5] := z[posb][1,2,3] * right[3,4,5]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef = coef1 * coef2
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	
	@tensor tmp2[6;4 1 5 2] := y[posb][1,2,3] * tmp1[3,4,5,6]	
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	
	@tensor threesitemps1[4,5;1,6,7,2] := x[posb][1,2,3] * tmp2[3,4,5,6,7]
	for (f1, f2) in fusiontrees(threesitemps1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 
		if coef != 1
			threesitemps1[f1, f2] .*= coef
		end			
	end

	# fuse indices
	cod = space(threesitemps1, 1) ⊗ space(threesitemps1, 2) 
	dom = space(threesitemps1, 3)' ⊗ space(threesitemps1, 4)'
	right4 = TensorMap(ds->zeros(scalartype(threesitemps1), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(threesitemps1)
		n = f2.uncoupled[2].n + f2.uncoupled[3].n + f2.uncoupled[4].n
		if n == 0
			f2′ = FusionTree((f2.uncoupled[1],f2.uncoupled[2]), f2.coupled, (f2.isdual[1],f2.isdual[2]))
			right4[f1, f2′] .+= threesitemps1[f1, f2][:,:,:,:,1,1]
		elseif n == 1
			isdual2 = ifelse(isodd(f2.uncoupled[2].n), f2.isdual[2], ifelse(isodd(f2.uncoupled[3].n), f2.isdual[3], f2.isdual[4] ))
			f2′ = FusionTree((f2.uncoupled[1], Z2Irrep(1)), f2.coupled, (f2.isdual[1], isdual2))
			right4[f1, f2′] .+= threesitemps1[f1, f2][:,:,:,:,1,1]
		end
	end

	@tensor tmp1[4;1 2 5 6] := z[posb-1][1,2,3] * right4[3,4,5,6]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef = coef1 * coef2
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	
	@tensor tmp2[6;4 1 5 2 7] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	
	@tensor threesitemps2[4,5;1,6,7,2,8] := x[posb-1][1,2,3] * tmp2[3,4,5,6,7,8]
	for (f1, f2) in fusiontrees(threesitemps2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 
		if coef != 1
			threesitemps2[f1, f2] .*= coef
		end			
	end		

	# fuse indices
	cod = space(threesitemps2, 1) ⊗ space(threesitemps2, 2) 
	dom = space(threesitemps2, 3)' ⊗ space(threesitemps2, 4)' ⊗ space(threesitemps2, 7)'
	right5 = TensorMap(ds->zeros(scalartype(threesitemps1), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(threesitemps2)
		n = f2.uncoupled[2].n + f2.uncoupled[3].n + f2.uncoupled[4].n
		if n == 0
			f2′ = FusionTree((f2.uncoupled[1],f2.uncoupled[2],f2.uncoupled[5]), f2.coupled, (f2.isdual[1],f2.isdual[2],f2.isdual[5]))
			right5[f1, f2′] .+= threesitemps2[f1, f2][:,:,:,:,1,1,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f2.uncoupled[2].n), f2.isdual[2], ifelse(isodd(f2.uncoupled[3].n), f2.isdual[3], f2.isdual[4] ))
			f2′ = FusionTree((f2.uncoupled[1], Z2Irrep(1), f2.uncoupled[5]), f2.coupled, (f2.isdual[1], isdual2, f2.isdual[5]))
			right5[f1, f2′] .+= threesitemps2[f1, f2][:,:,:,:,1,1,:]
		end
	end

	# trace physices
	right = TensorMap(ds->zeros(scalartype(right5), ds), space(right5, 1) ⊗ space(right5, 2) ← space(right5, 3)' )
	for (f1, f2) in fusiontrees(right5)
		if f2.uncoupled[2] == f2.uncoupled[3]
			f0 = FusionTree((f2.uncoupled[1],), f2.coupled, (f2.isdual[1],))
			@tensor right[f1, f0][1,2,3] += right5[f1, f2][1,2,3,4,4]
			# axpy!(1, tmp, r[f0, f2]) 
		end
	end	


	right = rmul!(right, f )

	return right
end
