# 5 GMPSs


function update_pair_left(left::AbstractTensorMap{<:ElementarySpace, 1, 5}, j::Int, x::Vector, y::Vector, z::Vector, u::Vector, v::Vector)
	posa = 2*j-1

	@tensor tmp1[1,2,3,4,7,8;5] := left[1,2,3,4,5,6] * v[posa][6,7,8] 
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end
	@tensor tmp2[1,2,3,5,8,6,9;4] := tmp1[1,2,3,4,5,6,7] * u[posa][7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end		
	tmp3 = g_fuse(tmp2, 4)

	@tensor tmp2[1,2,4,8,5,6,9;3] := tmp3[1,2,3,4,5,6,7] * z[posa][7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	tmp3 = g_fuse(tmp2, 3)	

	@tensor tmp2[1,3,8,4,5,6,9;2] := tmp3[1,2,3,4,5,6,7] * y[posa][7,8,9]	

	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	tmp3 = g_fuse(tmp2, 2)

	@tensor tmp2[1,2,8,3,4,5,6; 9] := tmp3[1,2,3,4,5,6,7] * x[posa][7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[3].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	
	tmp3 = g_fuse(tmp2, 2)

	@tensor tmp1[1,2,3,4,5,8,9;6] := tmp3[1,2,3,4,5,6,7] * x[posa+1][7,8,9]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end

	@tensor tmp2[1,2,3,4,6,9,7,10;5] := tmp1[1,2,3,4,5,6,7,8] * y[posa+1][8,9,10]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[8].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[7].n)) ? -1 : 1
		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end		
	tmp3 = g_fuse(tmp2, 5)

	@tensor tmp2[1,2,3,5,9,6,7,10;4] := tmp3[1,2,3,4,5,6,7,8] * z[posa+1][8,9,10]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[8].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[5].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end

	tmp3 = g_fuse(tmp2, 4)
	@tensor tmp2[1,2,4,9,5,6,7,10;3] := tmp3[1,2,3,4,5,6,7,8] * u[posa+1][8,9,10]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[6].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f1.uncoupled[7].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f1.uncoupled[8].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef8 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef9 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	tmp3 = g_fuse(tmp2, 3)
	@tensor tmp2[1,2,3,9,4,5,6,7;10] := tmp3[1,2,3,4,5,6,7,8] * v[posa+1][8,9,10]

	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f1.uncoupled[5].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[6].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[7].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1
		coef4 = (isodd(f1.uncoupled[8].n) && isodd(f1.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	tmp3 = g_fuse(tmp2, 3)

	left = g_trace(tmp3, 2)
	return permute(left, (1,), (2,3,4,5,6))
end

function update_pair_right(right::AbstractTensorMap{<:ElementarySpace, 5, 1}, j::Int, x::Vector, y::Vector, z::Vector, u::Vector, v::Vector)
	posb = 2 * j

	@tensor tmp1[4;5 6 7 1 2 8] := v[posb][1,2,3] * right[3,4,5,6,7,8]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[5].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef8 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	

	@tensor tmp2[4;5 6 1 7 2 8 9] := u[posb][1,2,3] * tmp1[3,4,5,6,7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[5].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end

	tmp3 = g_fuse(tmp2, 6)
	@tensor tmp2[4;5 1 6 7 2 8 9] := z[posb][1,2,3] * tmp3[3,4,5,6,7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[5].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	tmp3 = g_fuse(tmp2, 6)

	@tensor tmp2[4;1 5 6 7 2 8 9] := y[posb][1,2,3] * tmp3[3,4,5,6,7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[5].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	tmp3 = g_fuse(tmp2, 6)
	@tensor tmp2[1;4 5 6 7 2 8 9] := x[posb][1,2,3] * tmp3[3,4,5,6,7,8,9]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end

	tmp3 = g_fuse(tmp2, 6)
	@tensor tmp1[4; 1 2 5 6 7 8 9] := x[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9]
	for (f1, f2) in fusiontrees(tmp1)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1

		coef = coef1 * coef2 
		if coef != 1
			tmp1[f1, f2] .*= coef
		end
	end	

	@tensor tmp2[6;1 4 2 5 7 8 9 10] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7,8,9,10]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[3].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end
	tmp3 = g_fuse(tmp2, 4)

	@tensor tmp2[7;1 4 5 2 6 8 9 10] := z[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9,10]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[4].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[5].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	tmp3 = g_fuse(tmp2, 5)
	@tensor tmp2[8;1 4 5 6 2 7 9 10] := u[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9,10]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[5].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef5 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1
		coef6 = (isodd(f2.uncoupled[2].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef7 = (isodd(f2.uncoupled[3].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef8 = (isodd(f2.uncoupled[4].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1
		coef9 = (isodd(f2.uncoupled[6].n) && isodd(f1.uncoupled[1].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 * coef5 * coef6 * coef7 * coef8 * coef9
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	

	tmp3 = g_fuse(tmp2, 6)
	@tensor tmp2[1;4 5 6 7 2 8 9 10] := v[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9,10]
	for (f1, f2) in fusiontrees(tmp2)
		coef1 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef4 = (isodd(f2.uncoupled[5].n) && isodd(f2.uncoupled[4].n)) ? -1 : 1

		coef = coef1 * coef2 * coef3 * coef4 
		if coef != 1
			tmp2[f1, f2] .*= coef
		end
	end	
	tmp3 = g_fuse(tmp2, 6)

	right = g_trace(tmp3, 6)

	return permute(right, (1,2,3,4,5), (6,))
end



