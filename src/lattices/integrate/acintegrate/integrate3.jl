
function update_pair_left(left::GrassmannTensorMap{<:AbstractTensorMap{S, 1, 3}}, j::Int, x::Vector, y::Vector, z::Vector) where S
	posa = 2*j-1

	@tensor tmp1[6,4,5,3;2] := left[6,1,2,3] * x[posa][1,4,5] 
	@tensor tmp2[7,1,5,2,6;3] := tmp1[7,1,2,3,4] * y[posa][4,5,6]

	@tensor threesitemps1[8,1,2,6,3,4;7] := tmp2[8,1,2,3,4,5] * z[posa][5,6,7]
	threesitemps1 = get_data(threesitemps1)

	# fuse indices
	cod = space(threesitemps1, 1) ⊗ space(threesitemps1, 2) ⊗ space(threesitemps1, 5) ⊗ space(threesitemps1, 6)
	dom = space(threesitemps1, 7)'
	left4 = TensorMap(ds->zeros(scalartype(threesitemps1), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(threesitemps1)
		n = f1.uncoupled[2].n + f1.uncoupled[3].n + f1.uncoupled[4].n
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2],f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, (f1.isdual[1],f1.isdual[2],f1.isdual[5], f1.isdual[6]))
			left4[f1′, f2] .+= threesitemps1[f1, f2][:,:,1,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4] ))
			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6]), f1.coupled, (f1.isdual[1], isdual2, f1.isdual[5], f1.isdual[6]))
			left4[f1′, f2] .+= threesitemps1[f1, f2][:,:,1,1,:,:,:]
		end
	end
	left4 = GrassmannTensorMap(left4)

	@tensor tmp1[7,1,5,6,4;3] := left4[7,1,2,3,4] * x[posa+1][2,5,6]
	@tensor tmp2[8,1,2,6,3,7;4] := tmp1[8,1,2,3,4,5] * y[posa+1][5,6,7]	
	@tensor threesitemps2[9,1,2,3,7,4,5;8] := tmp2[9,1,2,3,4,5,6] * z[posa+1][6,7,8]


	# fuse indices
	threesitemps2 = get_data(threesitemps2)
	cod = space(threesitemps2, 1) ⊗ space(threesitemps2, 2) ⊗ space(threesitemps2, 3) ⊗ space(threesitemps2, 6) ⊗ space(threesitemps2, 7)
	dom = space(threesitemps2, 8)'
	left5 = TensorMap(ds->zeros(scalartype(threesitemps2), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(threesitemps2)
		n = f1.uncoupled[3].n + f1.uncoupled[4].n + f1.uncoupled[5].n
		if n == 0
			f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2],f1.uncoupled[3],f1.uncoupled[6], f1.uncoupled[7]), f1.coupled, 
							(f1.isdual[1],f1.isdual[2],f1.isdual[3],f1.isdual[6],f1.isdual[7]))
			left5[f1′, f2] .+= threesitemps2[f1, f2][:,:,:,1,1,:,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], ifelse(isodd(f1.uncoupled[4].n), f1.isdual[4], f1.isdual[5] ))
			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[6], f1.uncoupled[7]), f1.coupled, 
							(f1.isdual[1],f1.isdual[2], isdual2, f1.isdual[6], f1.isdual[7]))
			left5[f1′, f2] .+= threesitemps2[f1, f2][:,:,:,1,1,:,:,:]
		end
	end
	left5 = GrassmannTensorMap(left5)

	# trace physices
	return permute(g_trace(left5, 2), (1,), (2,3,4))
end

function update_pair_right(right::GrassmannTensorMap{<:AbstractTensorMap{S, 3, 1}}, j::Int, x::Vector, y::Vector, z::Vector) where S
	posb = 2 * j

	@tensor tmp1[4 ;1 2 5 6] := z[posb][1,2,3] * right[3,4,5,6]
	@tensor tmp2[6;4 1 5 2 7] := y[posb][1,2,3] * tmp1[3,4,5,6,7]	
	@tensor threesitemps1[4,5;1,6,7,2,8] := x[posb][1,2,3] * tmp2[3,4,5,6,7,8]

	# fuse indices
	threesitemps1 = get_data(threesitemps1)
	cod = space(threesitemps1, 1) ⊗ space(threesitemps1, 2) 
	dom = space(threesitemps1, 3)' ⊗ space(threesitemps1, 4)' ⊗ space(threesitemps1, 7)'
	right4 = TensorMap(ds->zeros(scalartype(threesitemps1), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(threesitemps1)
		n = f2.uncoupled[2].n + f2.uncoupled[3].n + f2.uncoupled[4].n
		if n == 0
			f2′ = FusionTree((f2.uncoupled[1],f2.uncoupled[2], f2.uncoupled[5]), f2.coupled, (f2.isdual[1],f2.isdual[2],f2.isdual[5]))
			right4[f1, f2′] .+= threesitemps1[f1, f2][:,:,:,:,1,1,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f2.uncoupled[2].n), f2.isdual[2], ifelse(isodd(f2.uncoupled[3].n), f2.isdual[3], f2.isdual[4] ))
			f2′ = FusionTree((f2.uncoupled[1], Z2Irrep(1),f2.uncoupled[5]), f2.coupled, (f2.isdual[1], isdual2, f2.isdual[5]))
			right4[f1, f2′] .+= threesitemps1[f1, f2][:,:,:,:,1,1,:]
		end
	end

	right4 = GrassmannTensorMap(right4)

	@tensor tmp1[4;1 2 5 6 7] := z[posb-1][1,2,3] * right4[3,4,5,6,7]
	@tensor tmp2[6;4 1 5 2 7 8] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7, 8]
	@tensor threesitemps2[4,5;1,6,7,2,8,9] := x[posb-1][1,2,3] * tmp2[3,4,5,6,7,8,9]

	# fuse indices
	threesitemps2 = get_data(threesitemps2)
	cod = space(threesitemps2, 1) ⊗ space(threesitemps2, 2) 
	dom = space(threesitemps2, 3)' ⊗ space(threesitemps2, 4)' ⊗ space(threesitemps2, 7)' ⊗ space(threesitemps2, 8)'
	right5 = TensorMap(ds->zeros(scalartype(threesitemps1), ds), cod ← dom) 
	for (f1, f2) in fusiontrees(threesitemps2)
		n = f2.uncoupled[2].n + f2.uncoupled[3].n + f2.uncoupled[4].n
		if n == 0
			f2′ = FusionTree((f2.uncoupled[1],f2.uncoupled[2],f2.uncoupled[5],f2.uncoupled[6]), f2.coupled, (f2.isdual[1],f2.isdual[2],f2.isdual[5],f2.isdual[6]))
			right5[f1, f2′] .+= threesitemps2[f1, f2][:,:,:,:,1,1,:,:]
		elseif n == 1
			isdual2 = ifelse(isodd(f2.uncoupled[2].n), f2.isdual[2], ifelse(isodd(f2.uncoupled[3].n), f2.isdual[3], f2.isdual[4] ))
			f2′ = FusionTree((f2.uncoupled[1], Z2Irrep(1), f2.uncoupled[5],f2.uncoupled[6]), f2.coupled, (f2.isdual[1], isdual2, f2.isdual[5],f2.isdual[6]))
			right5[f1, f2′] .+= threesitemps2[f1, f2][:,:,:,:,1,1,:,:]
		end
	end

	right5 = GrassmannTensorMap(right5)
	return permute(g_trace(right5, 4), (1,2,3), (4,))
end

