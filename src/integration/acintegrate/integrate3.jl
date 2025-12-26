
# function update_pair_left(left::GrassmannTensorMap{<:AbstractTensorMap{<:Number, S, 1, 3}}, j::Int, x::Vector, y::Vector, z::Vector) where S
# 	posa = 2*j-1

# 	@tensor tmp1[6,4,5,3;2] := left[6,1,2,3] * x[posa][1,4,5] 
# 	@tensor tmp2[7,1,5,2,6;3] := tmp1[7,1,2,3,4] * y[posa][4,5,6]

# 	@tensor threesitemps1[8,1,2,6,3,4;7] := tmp2[8,1,2,3,4,5] * z[posa][5,6,7]
# 	threesitemps1 = get_data(threesitemps1)

# 	# fuse indices
# 	cod = space(threesitemps1, 1) ⊗ space(threesitemps1, 2) ⊗ space(threesitemps1, 5) ⊗ space(threesitemps1, 6)
# 	dom = space(threesitemps1, 7)'
# 	left4 = zeros(scalartype(threesitemps1), cod ← dom) 
# 	for (f1, f2) in fusiontrees(threesitemps1)
# 		n = f1.uncoupled[2].n + f1.uncoupled[3].n + f1.uncoupled[4].n
# 		if n == 0
# 			f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2],f1.uncoupled[5], f1.uncoupled[6]), f1.coupled)
# 			left4[f1′, f2] .+= threesitemps1[f1, f2][:,:,1,1,:,:,:]
# 		elseif n == 1
# 			# isdual2 = ifelse(isodd(f1.uncoupled[2].n), f1.isdual[2], ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], f1.isdual[4] ))
# 			f1′ = FusionTree((f1.uncoupled[1], Z2Irrep(1), f1.uncoupled[5], f1.uncoupled[6]), f1.coupled)
# 			left4[f1′, f2] .+= threesitemps1[f1, f2][:,:,1,1,:,:,:]
# 		end
# 	end
# 	left4 = GrassmannTensorMap(left4)

# 	@tensor tmp1[7,1,5,6,4;3] := left4[7,1,2,3,4] * x[posa+1][2,5,6]
# 	@tensor tmp2[8,1,2,6,3,7;4] := tmp1[8,1,2,3,4,5] * y[posa+1][5,6,7]	
# 	@tensor threesitemps2[9,1,2,3,7,4,5;8] := tmp2[9,1,2,3,4,5,6] * z[posa+1][6,7,8]


# 	# fuse indices
# 	threesitemps2 = get_data(threesitemps2)
# 	cod = space(threesitemps2, 1) ⊗ space(threesitemps2, 2) ⊗ space(threesitemps2, 3) ⊗ space(threesitemps2, 6) ⊗ space(threesitemps2, 7)
# 	dom = space(threesitemps2, 8)'
# 	left5 = zeros(scalartype(threesitemps2), cod ← dom) 
# 	for (f1, f2) in fusiontrees(threesitemps2)
# 		n = f1.uncoupled[3].n + f1.uncoupled[4].n + f1.uncoupled[5].n
# 		if n == 0
# 			f1′ = FusionTree((f1.uncoupled[1], f1.uncoupled[2],f1.uncoupled[3],f1.uncoupled[6], f1.uncoupled[7]), f1.coupled)
# 			left5[f1′, f2] .+= threesitemps2[f1, f2][:,:,:,1,1,:,:,:]
# 		elseif n == 1
# 			# isdual2 = ifelse(isodd(f1.uncoupled[3].n), f1.isdual[3], ifelse(isodd(f1.uncoupled[4].n), f1.isdual[4], f1.isdual[5] ))
# 			f1′ = FusionTree((f1.uncoupled[1],f1.uncoupled[2], Z2Irrep(1), f1.uncoupled[6], f1.uncoupled[7]), f1.coupled)
# 			left5[f1′, f2] .+= threesitemps2[f1, f2][:,:,:,1,1,:,:,:]
# 		end
# 	end
# 	left5 = GrassmannTensorMap(left5)

# 	# trace physices
# 	return permute(g_trace(left5, 2), (1,), (2,3,4))
# end

# function update_pair_right(right::GrassmannTensorMap{<:AbstractTensorMap{<:Number, S, 3, 1}}, j::Int, x::Vector, y::Vector, z::Vector) where S
# 	posb = 2 * j

# 	@tensor tmp1[4 ;1 2 5 6] := z[posb][1,2,3] * right[3,4,5,6]
# 	@tensor tmp2[6;4 1 5 2 7] := y[posb][1,2,3] * tmp1[3,4,5,6,7]	
# 	@tensor threesitemps1[4,5;1,6,7,2,8] := x[posb][1,2,3] * tmp2[3,4,5,6,7,8]

# 	# fuse indices
# 	threesitemps1 = get_data(threesitemps1)
# 	cod = space(threesitemps1, 1) ⊗ space(threesitemps1, 2) 
# 	dom = space(threesitemps1, 3)' ⊗ space(threesitemps1, 4)' ⊗ space(threesitemps1, 7)'
# 	right4 = zeros(scalartype(threesitemps1), cod ← dom) 
# 	for (f1, f2) in fusiontrees(threesitemps1)
# 		n = f2.uncoupled[2].n + f2.uncoupled[3].n + f2.uncoupled[4].n
# 		if n == 0
# 			f2′ = FusionTree((f2.uncoupled[1],f2.uncoupled[2], f2.uncoupled[5]), f2.coupled)
# 			right4[f1, f2′] .+= threesitemps1[f1, f2][:,:,:,:,1,1,:]
# 		elseif n == 1
# 			# isdual2 = ifelse(isodd(f2.uncoupled[2].n), f2.isdual[2], ifelse(isodd(f2.uncoupled[3].n), f2.isdual[3], f2.isdual[4] ))
# 			f2′ = FusionTree((f2.uncoupled[1], Z2Irrep(1),f2.uncoupled[5]), f2.coupled)
# 			right4[f1, f2′] .+= threesitemps1[f1, f2][:,:,:,:,1,1,:]
# 		end
# 	end

# 	right4 = GrassmannTensorMap(right4)

# 	@tensor tmp1[4;1 2 5 6 7] := z[posb-1][1,2,3] * right4[3,4,5,6,7]
# 	@tensor tmp2[6;4 1 5 2 7 8] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7, 8]
# 	@tensor threesitemps2[4,5;1,6,7,2,8,9] := x[posb-1][1,2,3] * tmp2[3,4,5,6,7,8,9]

# 	# fuse indices
# 	threesitemps2 = get_data(threesitemps2)
# 	cod = space(threesitemps2, 1) ⊗ space(threesitemps2, 2) 
# 	dom = space(threesitemps2, 3)' ⊗ space(threesitemps2, 4)' ⊗ space(threesitemps2, 7)' ⊗ space(threesitemps2, 8)'
# 	right5 = zeros(scalartype(threesitemps1), cod ← dom) 
# 	for (f1, f2) in fusiontrees(threesitemps2)
# 		n = f2.uncoupled[2].n + f2.uncoupled[3].n + f2.uncoupled[4].n
# 		if n == 0
# 			f2′ = FusionTree((f2.uncoupled[1],f2.uncoupled[2],f2.uncoupled[5],f2.uncoupled[6]), f2.coupled)
# 			right5[f1, f2′] .+= threesitemps2[f1, f2][:,:,:,:,1,1,:,:]
# 		elseif n == 1
# 			# isdual2 = ifelse(isodd(f2.uncoupled[2].n), f2.isdual[2], ifelse(isodd(f2.uncoupled[3].n), f2.isdual[3], f2.isdual[4] ))
# 			f2′ = FusionTree((f2.uncoupled[1], Z2Irrep(1), f2.uncoupled[5],f2.uncoupled[6]), f2.coupled)
# 			right5[f1, f2′] .+= threesitemps2[f1, f2][:,:,:,:,1,1,:,:]
# 		end
# 	end

# 	right5 = GrassmannTensorMap(right5)
# 	return permute(g_trace(right5, 4), (1,2,3), (4,))
# end




function update_pair_left(left::GrassmannTensorMap{<:AbstractTensorMap{<:Number, S, 1, 3}}, j::Int, x::Vector, y::Vector, z::Vector) where S
	posa = 2*j-1

	@tensor tmp1[1,2,5,6; 3] := left[1,2,3,4] * z[posa][4,5,6] 
	@tensor tmp2[1,5,7,6,8;2] := tmp1[1,2,5,6,3] * y[posa][3,7,8]
	tmp3 = g_fuse(tmp2, 2)

	@tensor tmp4[1,57,9,6,8;a] := tmp3[1,57,6,8,2] * x[posa][2,9,a]
	tmp5 = g_fuse(tmp4, 2)

	# \bar{a}
	@tensor tmp1[1,579,f,g,6;8] := tmp5[1,579,6,8,a] * x[posa+1][a,f,g]
	@tensor tmp2[1,579,f,d,g,e;6] := tmp1[1,579,f,g,6,8] * y[posa+1][8,d,e]
	tmp3 = g_fuse(tmp2, 3)

	@tensor tmp4[1,579,fd,b;g,e,c] := tmp3[1,579,fd,g,e,6] * z[posa+1][6,b,c]
	tmp5 = g_fuse(tmp4, 3)	

	# trace physices
	left = g_trace(tmp5, 2)

	return left
end


function update_pair_right(right::GrassmannTensorMap{<:AbstractTensorMap{<:Number, S, 3, 1}}, j::Int, x::Vector, y::Vector, z::Vector) where S
	posb = 2 * j
	@tensor tmp1[2;9 a 3 4] := z[posb][9,a,1] * right[1,2,3,4]
	@tensor tmp2[3;9 7 a 8 4] := y[posb][7,8,2] * tmp1[2,9,a,3,4]
	tmp3 = g_fuse(tmp2, 4)

	@tensor tmp4[5;7 9 a8 6 4] := x[posb][5,6,3] * tmp3[3,9,7,a8,4]	
	tmp5 = g_fuse(tmp4, 4)

	@tensor tmp1[7;b c 9 a86 4] := x[posb-1][b,c,5] * tmp5[5,7,9,a86,4]
	@tensor tmp2[9;d b e c a86 4] := y[posb-1][d,e,7] * tmp1[7,b,c,9,a86,4]
	tmp3 = g_fuse(tmp2, 4)

	@tensor tmp4[f,d,b,g,ec,a86;4] := z[posb-1][f,g,9] * tmp3[9,d,b,ec,a86,4]	
	tmp5 = g_fuse(tmp4, 4)

	right = g_trace(tmp5, 4)	

	return right

end



# can be deleted
# function test()
# 	L = 8
# 	L2 = 2*L
# 	A = randomgmps(Float64, L2, D=6)
# 	B = randomgmps(Float64, L2, D=8)
# 	C = randomgmps(Float64, L2, D=10)
# 	m = GrassmannTransferMatrix(A,B,C)
# 	left = l_LL(randn, Z2Space(0=>1), m)
# 	right = r_RR(randn, Z2Space(0=>1), m)

# 	A = map(GrassmannTensorMap, A.data)
# 	B = map(GrassmannTensorMap, B.data)
# 	C = map(GrassmannTensorMap, C.data)
# 	for i in 1:L
# 		left1 = update_pair_left(left, i, A,B,C)
# 		left2 = update_pair_left2(left, i, A,B,C)
# 		println(norm(left1.data - left2.data))

# 		right1 = update_pair_right(right, L-i+1, A,B,C)
# 		right2 = update_pair_right2(right, L-i+1, A,B,C)
# 		println(norm(right1.data - right2.data))

# 		left = left1
# 		right = right1
# 	end
# end
# test();
