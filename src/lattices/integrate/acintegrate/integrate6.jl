

function update_pair_left(left::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS, u::GrassmannMPS, v::GrassmannMPS, w::GrassmannMPS; trunc)
	posa = 2*j-1
	f = (scaling(x) * scaling(y) * scaling(z) * scaling(u) * scaling(v) * scaling(w))^2

	@tensor tmp1[1,2,3,4,7,8;5] := left[1,2,3,4,5,6] * w[posa][6,7,8] 
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end
	@tensor tmp2[1,2,3,5,8,6,9;4] := tmp1[1,2,3,4,5,6,7] * v[posa][7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7)
	dom = space(tmp2, 8)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[4].n + f1.uncoupled[5].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[6], f1.uncoupled[7]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3], f1.isdual[4],f1.isdual[6],f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[4].n), f1.isdual[4], f1.isdual[5])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], Z2Irrep(1), f1.uncoupled[6], f1.uncoupled[7]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3], isdual2,f1.isdual[6], f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:]
		end
	end
	@tensor tmp2[1,2,4,8,5,6,9;3] := tmp3[1,2,3,4,5,6,7] * u[posa][7,8,9]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7)
	dom = space(tmp2, 8)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[3].n + f1.uncoupled[4].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[3], f1.isdual[5],f1.isdual[6],f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2],isdual2, f1.isdual[5], f1.isdual[6], f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:]
		end
	end
	@tensor tmp2[1,3,8,4,5,6,9;2] := tmp3[1,2,3,4,5,6,7] * z[posa][7,8,9]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7)
	dom = space(tmp2, 8)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[4], f1.isdual[5],f1.isdual[6],f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
							f1.coupled, (f1.isdual[1],isdual2,f1.isdual[4], f1.isdual[5], f1.isdual[6], f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:]
		end
	end
	@tensor tmp2[2,8,3,4,5,6,9;1] := tmp3[1,2,3,4,5,6,7] * y[posa][7,8,9]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef10 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef11 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) 
	dom = space(tmp2, 8)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[1].n + f1.uncoupled[2].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
								f1.coupled, (f1.isdual[1],f1.isdual[3],f1.isdual[4], f1.isdual[5],f1.isdual[6],f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[1].n), f1.isdual[1], f1.isdual[2])
			f1′ = FusionTree((Z2Irrep(1),f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
							f1.coupled, (isdual2,f1.isdual[3], f1.isdual[4], f1.isdual[5], f1.isdual[6], f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:,:,:]
		end
	end
	@tensor tmp2[1,8,2,3,4,5,6; 9] := tmp3[1,2,3,4,5,6,7] * x[posa][7,8,9]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end			
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7) 
	dom = space(tmp2, 8)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[1].n + f1.uncoupled[2].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
								f1.coupled, (f1.isdual[1],f1.isdual[3],f1.isdual[4], f1.isdual[5],f1.isdual[6],f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[1].n), f1.isdual[1], f1.isdual[2])
			f1′ = FusionTree((Z2Irrep(1),f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), 
							f1.coupled, (isdual2,f1.isdual[3], f1.isdual[4], f1.isdual[5], f1.isdual[6], f1.isdual[7]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:,:,:]
		end
	end

	# \bar{a}
	@tensor tmp1[1,2,3,4,5,8,9;6] := tmp3[1,2,3,4,5,6,7] * x[posa+1][7,8,9]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end

	@tensor tmp2[1,2,3,4,6,9,7,10;5] := tmp1[1,2,3,4,5,6,7,8] * y[posa+1][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5
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
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[4], f1.isdual[5],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,:,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[5].n), f1.isdual[5], f1.isdual[6])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], f1.uncoupled[4], Z2Irrep(1), f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[4], isdual2,f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,:,1,:,:,:]
		end
	end
	@tensor tmp2[1,2,3,5,9,6,7,10;4] := tmp3[1,2,3,4,5,6,7,8] * z[posa+1][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 
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
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[4], f1.isdual[6],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[4].n), f1.isdual[4], f1.isdual[5])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[3], Z2Irrep(1), f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], isdual2,f1.isdual[6], f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:,:]
		end
	end
	@tensor tmp2[1,2,4,9,5,6,7,10;3] := tmp3[1,2,3,4,5,6,7,8] * u[posa+1][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9
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
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[3], f1.isdual[5], f1.isdual[6],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4])
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (f1.isdual[1],f1.isdual[2],isdual2, f1.isdual[5], f1.isdual[6], f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:,:]
		end
	end
	@tensor tmp2[1,3,9,4,5,6,7,10;2] := tmp3[1,2,3,4,5,6,7,8] * v[posa+1][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[6].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[8].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef10 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef11 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9 * coef10 * coef11
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
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[4], f1.isdual[5], f1.isdual[6],f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
			f1′ = FusionTree((f1.uncoupled[1],Z2Irrep(1),f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7], f1.uncoupled[8]), 
							f1.coupled, (f1.isdual[1],isdual2, f1.isdual[4], f1.isdual[5], f1.isdual[6], f1.isdual[7], f1.isdual[8]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:]
		end
	end
	@tensor tmp2[1,2,9,3,4,5;6,7,10] := tmp3[1,2,3,4,5,6,7,8] * w[posa+1][8,9,10]
	# println("dim is ", dim(tmp2))
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end				
	# fuse physical indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6)  
	dom = space(tmp2, 7)' ⊗ space(tmp2, 8)' ⊗ space(tmp2, 9)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n 
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6]), 
								f1.coupled, (f1.isdual[1],f1.isdual[2], f1.isdual[4], f1.isdual[5], f1.isdual[6]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6]), 
							f1.coupled, (f1.isdual[1],isdual2, f1.isdual[4], f1.isdual[5], f1.isdual[6]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:,:,:]
		end
	end

	# trace physices
	left = TensorMap(ds->zeros(scalartype(tmp3), ds), space(tmp3, 3) ⊗ space(tmp3, 4) ⊗ space(tmp3, 5) ← space(tmp3, 6)' ⊗ space(tmp3, 7)' ⊗ space(tmp3, 8)')
	for (f1, f2) in fusiontrees(tmp3)
		if f1.uncoupled[1] == f1.uncoupled[2]
			f0 = FusionTree((f1.uncoupled[3],f1.uncoupled[4], f1.uncoupled[5]), 
				f1.coupled, (f1.isdual[3],f1.isdual[4],f1.isdual[5]))
			@tensor left[f0, f2][2,3,4,5,6,7] += tmp3[f1, f2][1,1,2,3,4,5,6,7]
			# axpy!(1, tmp, r[f0, f2]) 
		end
	end	

	left = rmul!(left, f)	

	return left	
end
