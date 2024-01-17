

function update_pair_left(left::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS, 
							z::GrassmannMPS, u::GrassmannMPS, v::GrassmannMPS, w::GrassmannMPS, q::GrassmannMPS; trunc=DefaultIntegrationTruncation)
	posa = 2*j-1
	f = (scaling(x) * scaling(y) * scaling(z) * scaling(u) * scaling(v) * scaling(w) * scaling(q))^2

	@tensor tmp1[1,2,3,4,5,8,9;6] := left[1,2,3,4,5,6,7] * q[posa][7,8,9] 
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		# println(coef1, " ", coef2, " ", coef3, " ", coef4, " ", coef5)
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end
	@tensor tmp2[1,2,3,4,9,6,10,7;5] := tmp1[1,2,3,4,5,6,7,8] * w[posa][8,9,10]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8)
	dom = space(tmp2, 9)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[5].n + f1.uncoupled[6].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[7], f1.uncoupled[8]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3], f1.isdual[4],f1.isdual[5],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,:,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[5].n), f1.isdual[5], f1.isdual[6])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], Z2Irrep(1), f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3], f1.isdual[4],isdual2,f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,:,1,:,:,:]
		end
	end
	@tensor tmp2[1,2,3,9,5,10,6,7;4] := tmp3[1,2,3,4,5,6,7,8] * v[posa][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef10 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8)
	dom = space(tmp2, 9)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[4].n + f1.uncoupled[5].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3], f1.isdual[4],f1.isdual[6],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[4].n), f1.isdual[4], f1.isdual[5])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], Z2Irrep(1), f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3],isdual2, f1.isdual[6], f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:,:]
		end
	end
	@tensor tmp2[1,2,9,4,10,5,6,7;3] := tmp3[1,2,3,4,5,6,7,8] * u[posa][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef10 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef11 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1			
		coef12 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef13 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1			
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11 * coef12 * coef13
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8)
	dom = space(tmp2, 9)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[3].n + f1.uncoupled[4].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3], f1.isdual[5],f1.isdual[6],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2], isdual2, f1.isdual[5], f1.isdual[6], f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:,:]
		end
	end
	@tensor tmp2[1,9,3,10,4,5,6,7;2] := tmp3[1,2,3,4,5,6,7,8] * z[posa][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef10 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef11 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef12 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef13 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef14 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef15 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef16 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11 * coef12 * coef13 * coef14 * coef15 * coef16
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8)
	dom = space(tmp2, 9)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[4], f1.isdual[5],f1.isdual[6],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (f1.isdual[1],isdual2, f1.isdual[4], f1.isdual[5], f1.isdual[6], f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:]
		end
	end
	@tensor tmp2[9,2,10,3,4,5,6,7; 1] := tmp3[1,2,3,4,5,6,7,8] * y[posa][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef8 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef10 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef11 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef12 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef13 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef14 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef15 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef16 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef17 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef18 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef19 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11 * coef12 * coef13 * coef14 * coef15 * coef16 * coef17 * coef18 * coef19
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8)
	dom = space(tmp2, 9)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[1].n + f1.uncoupled[2].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
								f1.coupled, (f1.isdual[1],f1.isdual[3],f1.isdual[4], f1.isdual[5],f1.isdual[6],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[1].n), f1.isdual[1], f1.isdual[2])
			f1′ = FusionTree((Z2Irrep(1), f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (isdual2, f1.isdual[3], f1.isdual[4], f1.isdual[5], f1.isdual[6], f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:,:,:,:]
		end
	end
	@tensor tmp2[9,1,10,2,3,4;5,6,7] := tmp3[1,2,3,4,5,6,7,8] * x[posa][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1			
		coef6 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1			
		coef8 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef9 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1			
		coef10 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef11 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1			
		coef12 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef13 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1			
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11 * coef12 * coef13 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) 
	dom = space(tmp2, 7)' ⊗ space(tmp2, 8)' ⊗ space(tmp2, 9)'
	mpsj1 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[1].n + f1.uncoupled[2].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6]), 
								f1.coupled, (f1.isdual[1],f1.isdual[3],f1.isdual[4], f1.isdual[5],f1.isdual[6]))
			mpsj1[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[1].n), f1.isdual[1], f1.isdual[2])
			f1′ = FusionTree((Z2Irrep(1), f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6]), 
							f1.coupled, (isdual2, f1.isdual[3], f1.isdual[4], f1.isdual[5], f1.isdual[6]))
			mpsj1[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:,:,:,:]
		end
	end

	# \bar{a}

	@tensor tmp1[1,2,3,4,5,6,9,10;7] := mpsj1[1,2,3,4,5,6,7,8] * q[posa+1][8,9,10]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end

	@tensor tmp2[1,2,3,4,5,7,10,8,11;6] := tmp1[1,2,3,4,5,6,7,8,9] * w[posa+1][9,10,11]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[9].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end		
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 8) ⊗ space(tmp2, 9) 
	dom = space(tmp2, 10)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[6].n + f1.uncoupled[7].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[8], f1.uncoupled[9]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[4], f1.isdual[5],f1.isdual[6], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,:,:,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[6].n), f1.isdual[6], f1.isdual[7])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], Z2Irrep(1), f1.uncoupled[8], f1.uncoupled[9]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[4], f1.isdual[5], isdual2, f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,:,:,1,:,:,:]
		end
	end
	@tensor tmp2[1,2,3,4,6,10,7,8,11;5] := tmp3[1,2,3,4,5,6,7,8,9] * v[posa+1][9,10,11]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[9].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end				
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8) ⊗ space(tmp2, 9) 
	dom = space(tmp2, 10)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[5].n + f1.uncoupled[6].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[7], f1.uncoupled[8], f1.uncoupled[9]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[4], f1.isdual[5],f1.isdual[7], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,:,1,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[5].n), f1.isdual[5], f1.isdual[6])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], Z2Irrep(1), f1.uncoupled[7], f1.uncoupled[8], f1.uncoupled[9]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[4], isdual2, f1.isdual[7], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,:,1,:,:,:,:]
		end
	end
	@tensor tmp2[1,2,3,5,10,6,7,8,11;4] := tmp3[1,2,3,4,5,6,7,8,9] * u[posa+1][9,10,11]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[9].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end				
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8) ⊗ space(tmp2, 9) 
	dom = space(tmp2, 10)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[4].n + f1.uncoupled[5].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8], f1.uncoupled[9]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[4], f1.isdual[6],f1.isdual[7], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[4].n), f1.isdual[4], f1.isdual[5])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], Z2Irrep(1), f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8], f1.uncoupled[9]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], isdual2, f1.isdual[6], f1.isdual[7], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:,:,:]
		end
	end
	@tensor tmp2[1,2,4,10,5,6,7,8,11;3] := tmp3[1,2,3,4,5,6,7,8,9] * z[posa+1][9,10,11]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[9].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef10 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef11 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end				
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8) ⊗ space(tmp2, 9) 
	dom = space(tmp2, 10)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[3].n + f1.uncoupled[4].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8], f1.uncoupled[9]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[5], f1.isdual[6],f1.isdual[7], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8], f1.uncoupled[9]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2], isdual2, f1.isdual[5], f1.isdual[6], f1.isdual[7], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:,:,:]
		end
	end
	@tensor tmp2[1,3,10,4,5,6,7,8,11;2] := tmp3[1,2,3,4,5,6,7,8,9] * y[posa+1][9,10,11]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[9].n)) ? -1 : 1
		coef8 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef10 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef11 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef12 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef13 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11 * coef12 * coef13
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end				
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) ⊗ space(tmp2, 8) ⊗ space(tmp2, 9) 
	dom = space(tmp2, 10)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8], f1.uncoupled[9]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[4], f1.isdual[5], f1.isdual[6],f1.isdual[7], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8], f1.uncoupled[9]), 
							f1.coupled, (f1.isdual[1],isdual2, f1.isdual[4], f1.isdual[5], f1.isdual[6], f1.isdual[7], f1.isdual[8],f1.isdual[9]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:,:]
		end
	end
	@tensor tmp2[1,2,10,11,8,7,6;5,4,3] := tmp3[1,2,3,4,5,6,7,8,9] * x[posa+1][9,10,11]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1

		coef8 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef9 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef10 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef11 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef12 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef13 = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1			

		coef14 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef15 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef16 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef17 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef18 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1

		coef19 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef20 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef21 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef22 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1

		coef23 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef24 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef25 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1

		coef26 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef27 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11 * coef12 * coef13 * coef14
		coef = coef * coef15 * coef16 * coef17 * coef18 * coef19 * coef20 * coef21 * coef22 * coef23 * coef24 * coef25 * coef26 * coef27
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end						

	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) 
	dom = space(tmp2, 8)' ⊗ space(tmp2, 9)' ⊗ space(tmp2, 10)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[4], f1.isdual[5], f1.isdual[6],f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
							f1.coupled, (f1.isdual[1],isdual2, f1.isdual[4], f1.isdual[5], f1.isdual[6], f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:,:]
		end
	end

	# trace physices
	left = TensorMap(ds->zeros(scalartype(tmp3), ds), space(tmp3, 3) ⊗ space(tmp3, 4) ⊗ space(tmp3, 5) ⊗ space(tmp3, 6) ← space(tmp3, 7)' ⊗ space(tmp3, 8)' ⊗ space(tmp3, 9)')
	for (f1, f2) in fusiontrees(tmp3)
		if f1.uncoupled[1] == f1.uncoupled[2]
			f0 = FusionTree((f1.uncoupled[3],f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6]), 
				f1.coupled, (f1.isdual[3],f1.isdual[4],f1.isdual[5],f1.isdual[6]))
			@tensor left[f0, f2][2,3,4,5,6,7,8] += tmp3[f1, f2][1,1,2,3,4,5,6,7,8]
			# axpy!(1, tmp, r[f0, f2]) 
		end
	end	

	left = rmul!(left, f)	

	return left	
end
