# 1 gmps integration


# function update_pair_left(left::AbstractTensorMap, j::Int, x::GrassmannMPS; trunc=DefaultIntegrationTruncation)
# 	f = scaling(x)^2
# 	pos1, pos2 = 2*j-1, 2*j
# 	@tensor tmp2[2,4; 5] := left[1] * x[pos1][1,2,3] * x[pos2][3,4,5]
# 	# fuse indices
# 	dom = space(tmp2, 3)'
# 	cod = one(dom)
# 	tmp1 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
# 	for (f1, f2) in fusiontrees(tmp2)
# 		if f1.uncoupled[1] == f1.uncoupled[2]
# 			f0 = FusionTree((), Z2Irrep(0), ())
# 			@tensor tmp1[f0, f2][2] += tmp2[f1, f2][1,1,2]
# 		end
# 	end	
# 	return rmul!(tmp1, f)
# end


# function update_pair_right(right::AbstractTensorMap, j::Int, x::GrassmannMPS; trunc=DefaultIntegrationTruncation)
# 	f = scaling(x)^2
# 	pos1, pos2 = 2*j-1, 2*j
# 	@tensor tmp2[1; 2 4] := x[pos1][1,2,3] * (x[pos2][3,4,5] * right[5])
# 	# fuse indices
# 	cod = space(tmp2, 1)
# 	dom = one(cod)
# 	tmp1 = TensorMap(ds->zeros(scalartype(tmp2), ds), cod ← dom) 
# 	for (f1, f2) in fusiontrees(tmp2)
# 		if f2.uncoupled[1] == f2.uncoupled[2]
# 			f0 = FusionTree((), Z2Irrep(0), ())
# 			@tensor tmp1[f1, f0][2] += tmp2[f1, f2][2,1,1]
# 		end
# 	end	
# 	return rmul!(tmp1, f)
# end

function update_pair_left(left::MPSBondTensor, j::Int, x::Vector; trunc=DefaultIntegrationTruncation)
	pos1, pos2 = 2*j-1, 2*j
	@tensor tmp2[-1, -2,-3; -4] := left[-1, 1] * x[pos1][1,-2,3] * x[pos2][3,-3,-4]
	return _trace_physical(tmp2)
end

function update_pair_right(right::MPSBondTensor, j::Int, x::Vector; trunc=DefaultIntegrationTruncation)
	pos1, pos2 = 2*j-1, 2*j
	@tensor tmp2[-1, -2,-3; -4] := x[pos1][-1,-2,1] * x[pos2][1,-3,2] * right[2,-4]
	return _trace_physical(tmp2)
end