
# function update_pair_left(left::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS, z::GrassmannMPS, u::GrassmannMPS; trunc=DefaultIntegrationTruncation)
# 	posa = 2*j-1
# 	f = (scaling(x) * scaling(y) * scaling(z) * scaling(u))^2	

# 	@tensor tmp1[1,2,5,6; 3] := left[1,2,3,4] * u[posa][4,5,6] 
# 	for (f1, f2) in fusiontrees(tmp1)
# 		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef2 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		# println(coef1, " ", coef2, " ", coef3, " ", coef4, " ", coef5)
# 		coef = coef1 * coef2 
# 		if coef != 1
# 			tmp1[f1, f2] .*= coef
# 		end
# 	end
# 	@tensor tmp2[1,3,6,4,7;2] := tmp1[1,2,3,4,5] * z[posa][5,6,7]
# 	for (f1, f2) in fusiontrees(tmp2)
# 		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef2 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef3 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef4 = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef5 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
# 		coef = coef1 * coef2 * coef3 * coef4 * coef5
# 		if coef != 1
# 			tmp2[f1, f2] .*= coef
# 		end
# 	end		

# 	# fuse indices
# 	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5)
# 	dom = space(tmp2, 6)'
# 	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
# 	for (f1, f2) in fusiontrees(tmp2)
# 		n = f1.uncoupled[2].n + f1.uncoupled[3].n
# 		if n == 0
# 			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[4],f1.isdual[5]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:]
# 		elseif n == 1
# 			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
# 			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, (f1.isdual[1], isdual2, f1.isdual[4], f1.isdual[5]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:]
# 		end
# 	end

# 	@tensor tmp2[2,6,3,4,7; 1] := tmp3[1,2,3,4,5] * y[posa][5,6,7]
# 	for (f1, f2) in fusiontrees(tmp2)
# 		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef2 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef3 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef4 = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef5 = (isodd(f1.uncoupled[1].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef6 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
# 		coef7 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
# 		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
# 		if coef != 1
# 			tmp2[f1, f2] .*= coef
# 		end
# 	end

# 	# fuse indices
# 	cod = space(tmp2, 1) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5)
# 	dom = space(tmp2, 6)'
# 	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
# 	for (f1, f2) in fusiontrees(tmp2)
# 		n = f1.uncoupled[1].n + f1.uncoupled[2].n
# 		if n == 0
# 			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[3],f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, (f1.isdual[1],f1.isdual[3],f1.isdual[4],f1.isdual[5]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:]
# 		elseif n == 1
# 			isdual2 = ifelse(isodd(f1.uncoupled[1].n), f1.isdual[1], f1.isdual[2])
# 			f1′ = FusionTree((Z2Irrep(1), f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, (isdual2, f1.isdual[3], f1.isdual[4], f1.isdual[5]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:]
# 		end
# 	end

# 	@tensor tmp2[1,6,2,3,4;7] := tmp3[1,2,3,4,5] * x[posa][5,6,7]
# 	for (f1, f2) in fusiontrees(tmp2)
# 		coef1 = (isodd(f1.uncoupled[3].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
# 		coef2 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
# 		coef3 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1
# 		coef = coef1 * coef2 * coef3 
# 		if coef != 1
# 			tmp2[f1, f2] .*= coef
# 		end
# 	end

# 	# fuse indices
# 	cod = space(tmp2, 1) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5)
# 	dom = space(tmp2, 6)'
# 	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
# 	for (f1, f2) in fusiontrees(tmp2)
# 		n = f1.uncoupled[1].n + f1.uncoupled[2].n
# 		if n == 0
# 			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[3],f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, (f1.isdual[1],f1.isdual[3],f1.isdual[4],f1.isdual[5]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:]
# 		elseif n == 1
# 			isdual2 = ifelse(isodd(f1.uncoupled[1].n), f1.isdual[1], f1.isdual[2])
# 			f1′ = FusionTree((Z2Irrep(1), f1.uncoupled[3], f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, (isdual2, f1.isdual[3], f1.isdual[4], f1.isdual[5]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,1,:,:,:,:]
# 		end
# 	end

# 	# \bar{a}
# 	@tensor tmp1[1,2,3,6,7;4] := tmp3[1,2,3,4,5] * x[posa+1][5,6,7]
# 	for (f1, f2) in fusiontrees(tmp1)
# 		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef = coef1 * coef2 
# 		if coef != 1
# 			tmp1[f1, f2] .*= coef
# 		end
# 	end	

# 	@tensor tmp2[1,2,4,7,5,8;3] := tmp1[1,2,3,4,5,6] * y[posa+1][6,7,8]
# 	for (f1, f2) in fusiontrees(tmp2)
# 		coef1 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef3 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef4 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef5 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
# 		coef = coef1 * coef2 * coef3 * coef4 * coef5
# 		if coef != 1
# 			tmp2[f1, f2] .*= coef
# 		end
# 	end	

# 	# fuse indices
# 	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6)
# 	dom = space(tmp2, 7)'
# 	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
# 	for (f1, f2) in fusiontrees(tmp2)
# 		n = f1.uncoupled[3].n + f1.uncoupled[4].n
# 		if n == 0
# 			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[3],f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, 
# 				(f1.isdual[1],f1.isdual[2],f1.isdual[3],f1.isdual[5],f1.isdual[6]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:]
# 		elseif n == 1
# 			isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4])
# 			f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, 
# 				(f1.isdual[1], f1.isdual[2], isdual2, f1.isdual[5], f1.isdual[6]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:]
# 		end
# 	end		

# 	@tensor tmp2[1,3,7,4,5,8;2] := tmp3[1,2,3,4,5,6] * z[posa+1][6,7,8]
# 	for (f1, f2) in fusiontrees(tmp2)
# 		coef1 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef3 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef4 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef5 = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
# 		coef6 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
# 		coef7 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
# 		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
# 		if coef != 1
# 			tmp2[f1, f2] .*= coef
# 		end
# 	end	

# 	# fuse indices
# 	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6)
# 	dom = space(tmp2, 7)'
# 	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
# 	for (f1, f2) in fusiontrees(tmp2)
# 		n = f1.uncoupled[2].n + f1.uncoupled[3].n
# 		if n == 0
# 			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[4],f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, 
# 				(f1.isdual[1],f1.isdual[2],f1.isdual[4],f1.isdual[5],f1.isdual[6]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:]
# 		elseif n == 1
# 			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
# 			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, 
# 				(f1.isdual[1], isdual2, f1.isdual[4], f1.isdual[5], f1.isdual[6]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:]
# 		end
# 	end		

# 	@tensor tmp2[1,2,7,3,4;5,8] := tmp3[1,2,3,4,5,6] * u[posa+1][6,7,8]
# 	for (f1, f2) in fusiontrees(tmp2)
# 		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
# 		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
# 		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
# 		coef = coef1 * coef2 * coef3
# 		if coef != 1
# 			tmp2[f1, f2] .*= coef
# 		end
# 	end	

# 	# fuse indices
# 	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) 
# 	dom = space(tmp2, 6)' ⊗ space(tmp2, 7)'
# 	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
# 	for (f1, f2) in fusiontrees(tmp2)
# 		n = f1.uncoupled[2].n + f1.uncoupled[3].n
# 		if n == 0
# 			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[4],f1.uncoupled[5]), f1.coupled, 
# 				(f1.isdual[1],f1.isdual[2],f1.isdual[4],f1.isdual[5]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:]
# 		elseif n == 1
# 			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
# 			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5]), f1.coupled, 
# 				(f1.isdual[1], isdual2, f1.isdual[4], f1.isdual[5]))
# 			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:]
# 		end
# 	end		

# 	# trace physices
# 	left = TensorMap(ds->zeros(scalartype(tmp3), ds), space(tmp3, 3) ⊗ space(tmp3, 4) ← space(tmp3, 5)' ⊗ space(tmp3, 6)')
# 	for (f1, f2) in fusiontrees(tmp3)
# 		if f1.uncoupled[1] == f1.uncoupled[2]
# 			f0 = FusionTree((f1.uncoupled[3],f1.uncoupled[4]), f1.coupled, (f1.isdual[3],f1.isdual[4]))
# 			@tensor left[f0, f2][2,3,4,5] += tmp3[f1, f2][1,1,2,3,4,5]
# 			# axpy!(1, tmp, r[f0, f2]) 
# 		end
# 	end	

# 	left = rmul!(left, f)	

# 	return left	
# end



function update_pair_left(left::AbstractTensorMap{<:ElementarySpace, 1, 4}, j::Int, x::Vector, y::Vector, z::Vector, u::Vector)
	posa = 2*j-1

	@tensor tmp1[7,1,2,5,6; 3] := left[7,1,2,3,4] * u[posa][4,5,6] 
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		# println(coef1, " ", coef2, " ", coef3, " ", coef4, " ", coef5)
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end
	@tensor tmp2[8,1,3,6,4,7;2] := tmp1[8,1,2,3,4,5] * z[posa][5,6,7]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end		

	# fuse indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6)
	dom = space(tmp2, 7)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[3].n + f1.uncoupled[4].n
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[3], f1.uncoupled[5],f1.uncoupled[6]), f1.coupled, 
							(f1.isdual[1],f1.isdual[2],f1.isdual[3],f1.isdual[5], f1.isdual[6]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4])
			f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, 
						(f1.isdual[1], f1.isdual[2], isdual2, f1.isdual[5], f1.isdual[6]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:]
		end
	end

	@tensor tmp2[8,2,6,3,4,7; 1] := tmp3[8,1,2,3,4,5] * y[posa][5,6,7]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end

	# fuse indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6)
	dom = space(tmp2, 7)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[4], f1.uncoupled[5],f1.uncoupled[6]), f1.coupled, 
							(f1.isdual[1],f1.isdual[2],f1.isdual[4],f1.isdual[5], f1.isdual[6]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, 
							(f1.isdual[1], isdual2, f1.isdual[4], f1.isdual[5], f1.isdual[6]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:]
		end
	end

	@tensor tmp2[8,1,6,2,3,4;7] := tmp3[8,1,2,3,4,5] * x[posa][5,6,7]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end

	# fuse indices
	cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 4) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6)
	dom = space(tmp2, 7)'
	tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(tmp2)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[4], f1.uncoupled[5],f1.uncoupled[6]), f1.coupled, 
							(f1.isdual[1],f1.isdual[2],f1.isdual[4],f1.isdual[5],f1.isdual[6]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], f1.isdual[3])
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[4], f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, 
							(f1.isdual[1], isdual2, f1.isdual[4], f1.isdual[5], f1.isdual[6]))
			tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,1,:,:,:,:]
		end
	end

	# \bar{a}
	@tensor tmp1[8,1,2,3,6,7;4] := tmp3[8,1,2,3,4,5] * x[posa+1][5,6,7]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	

	@tensor tmp2[9,1,2,4,7,5,8;3] := tmp1[9,1,2,3,4,5,6] * y[posa+1][6,7,8]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	# fuse indices
	# cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 4) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7)
	# dom = space(tmp2, 8)'
	# tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	# for (f1, f2) in fusiontrees(tmp2)
	# 	n = f1.uncoupled[4].n + f1.uncoupled[5].n
	# 	if n == 0
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[3],f1.uncoupled[4], f1.uncoupled[6],f1.uncoupled[7]), f1.coupled, 
	# 			(f1.isdual[1],f1.isdual[2],f1.isdual[3],f1.isdual[4],f1.isdual[6],f1.isdual[7]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:]
	# 	elseif n == 1
	# 		isdual2 = ifelse(isodd(f1.uncoupled[4].n), f1.isdual[4], f1.isdual[5])
	# 		f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2],f1.uncoupled[3], Z2Irrep(1), f1.uncoupled[6], f1.uncoupled[7]), f1.coupled, 
	# 			(f1.isdual[1], f1.isdual[2], f1.isdual[3], isdual2, f1.isdual[6], f1.isdual[7]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,:,1,:,:,:]
	# 	end
	# end		
	tmp3 = g_fuse(tmp2, 4)

	@tensor tmp2[9,1,3,7,4,5,8;2] := tmp3[9,1,2,3,4,5,6] * z[posa+1][6,7,8]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	# fuse indices
	# cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) ⊗ space(tmp2, 7)
	# dom = space(tmp2, 8)'
	# tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	# for (f1, f2) in fusiontrees(tmp2)
	# 	n = f1.uncoupled[3].n + f1.uncoupled[4].n
	# 	if n == 0
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[3],f1.uncoupled[5], f1.uncoupled[6],f1.uncoupled[7]), f1.coupled, 
	# 			(f1.isdual[1],f1.isdual[2],f1.isdual[3],f1.isdual[5],f1.isdual[6],f1.isdual[7]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:]
	# 	elseif n == 1
	# 		isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4])
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6], f1.uncoupled[7]), f1.coupled, 
	# 			(f1.isdual[1],f1.isdual[2], isdual2, f1.isdual[5], f1.isdual[6], f1.isdual[7]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:]
	# 	end
	# end	
	tmp3 = g_fuse(tmp2, 3)	

	@tensor tmp2[9,1,2,7,3,4;5,8] := tmp3[9,1,2,3,4,5,6] * u[posa+1][6,7,8]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	# fuse indices
	# cod = space(tmp2, 1) ⊗ space(tmp2, 2) ⊗ space(tmp2, 3) ⊗ space(tmp2, 5) ⊗ space(tmp2, 6) 
	# dom = space(tmp2, 7)' ⊗ space(tmp2, 8)'
	# tmp3 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
	# for (f1, f2) in fusiontrees(tmp2)
	# 	n = f1.uncoupled[3].n + f1.uncoupled[4].n
	# 	if n == 0
	# 		f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2],f1.uncoupled[3],f1.uncoupled[5],f1.uncoupled[6]), f1.coupled, 
	# 			(f1.isdual[1],f1.isdual[2],f1.isdual[3],f1.isdual[5],f1.isdual[6]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:]
	# 	elseif n == 1
	# 		isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4])
	# 		f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, 
	# 			(f1.isdual[1], f1.isdual[2], isdual2, f1.isdual[5], f1.isdual[6]))
	# 		tmp3[f1′, f2] .+= tmp2[f1, f2][:,:,:,1,:,:,:,:]
	# 	end
	# end		
	tmp3 = g_fuse(tmp2, 3)

	# trace physices
	# left = TensorMap(ds->zeros(scalartype(tmp3), ds), space(tmp3, 1) ⊗ space(tmp3, 4) ⊗ space(tmp3, 5) ← space(tmp3, 6)' ⊗ space(tmp3, 7)')
	# for (f1, f2) in fusiontrees(tmp3)
	# 	if f1.uncoupled[2] == f1.uncoupled[3]
	# 		f0 = FusionTree((f1.uncoupled[1],f1.uncoupled[4],f1.uncoupled[5]), f1.coupled, (f1.isdual[1],f1.isdual[4],f1.isdual[5]))
	# 		@tensor left[f0, f2][1,3,4,5,6] += tmp3[f1, f2][1,2,2,3,4,5,6]
	# 		# axpy!(1, tmp, r[f0, f2]) 
	# 	end
	# end	
	left = g_trace(tmp3, 2)

	return permute(left, (1,), (2,3,4,5))
end


function update_pair_right(right::AbstractTensorMap{<:ElementarySpace, 4, 1}, j::Int, x::Vector, y::Vector, z::Vector, u::Vector)
	posb = 2 * j
	@tensor tmp1[4 ;5 6 1 2 7] := u[posb][1,2,3] * right[3,4,5,6,7]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	
	@tensor tmp2[4;5 1 6 2 7 8] := z[posb][1,2,3] * tmp1[3,4,5,6,7,8]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end

	# fuse indices
	tmp3 = g_fuse(tmp2, 5)
	# println(space(y[posb], 3), " ", space(tmp2,1), " ", space(tmp3,1), " ", space(right, 3))

	@tensor tmp2[4;1 5 6 2 7 8] := y[posb][1,2,3] * tmp3[3,4,5,6,7,8]	
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	# fuse indices
	tmp3 = g_fuse(tmp2, 5)

	@tensor tmp2[1 ;4 5 6 2 7 8] := x[posb][1,2,3] * tmp3[3,4,5,6,7,8]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end		

	# fuse indices
	tmp3 = g_fuse(tmp2, 5)

	@tensor tmp1[4 ;5 6 1 2 7 8] := x[posb-1][1,2,3] * tmp3[3,4,5,6,7,8]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	

	@tensor tmp2[4;5 1 6 2 7 8 9] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	tmp3 = g_fuse(tmp2, 5)

	@tensor tmp2[4;1 5 6 2 7 8 9] := z[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9]	
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	
	tmp3 = g_fuse(tmp2, 5)

	@tensor tmp2[1 ;4 5 6 2 7 8 9] := u[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	
	tmp3 = g_fuse(tmp2, 5)

	right = g_trace(tmp3, 5)	

	return permute(right, (1,2,3,4), (5,))

end
