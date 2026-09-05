# 5 GMPSs


function update_pair_left(left::AbstractParityTensorMap{<:Number, 1, 5}, j::Int, x::Vector, y::Vector, z::Vector, u::Vector, v::Vector)
	posa = 2*j-1

	@grassmann tmp1[1,2,3,4,7,8;5] := left[1,2,3,4,5,6] * v[posa][6,7,8]
	@grassmann tmp2[1,2,3,5,8,6,9;4] := tmp1[1,2,3,4,5,6,7] * u[posa][7,8,9]
	tmp3 = g_fuse(tmp2, 4)

	@grassmann tmp2[1,2,4,8,5,6,9;3] := tmp3[1,2,3,4,5,6,7] * z[posa][7,8,9]
	tmp3 = g_fuse(tmp2, 3)

	@grassmann tmp2[1,3,8,4,5,6,9;2] := tmp3[1,2,3,4,5,6,7] * y[posa][7,8,9]
	tmp3 = g_fuse(tmp2, 2)

	@grassmann tmp2[1,2,8,3,4,5,6; 9] := tmp3[1,2,3,4,5,6,7] * x[posa][7,8,9]
	tmp3 = g_fuse(tmp2, 2)

	@grassmann tmp1[1,2,3,4,5,8,9;6] := tmp3[1,2,3,4,5,6,7] * x[posa+1][7,8,9]
	@grassmann tmp2[1,2,3,4,6,9,7,10;5] := tmp1[1,2,3,4,5,6,7,8] * y[posa+1][8,9,10]

	tmp3 = g_fuse(tmp2, 5)

	@grassmann tmp2[1,2,3,5,9,6,7,10;4] := tmp3[1,2,3,4,5,6,7,8] * z[posa+1][8,9,10]

	tmp3 = g_fuse(tmp2, 4)
	@grassmann tmp2[1,2,4,9,5,6,7,10;3] := tmp3[1,2,3,4,5,6,7,8] * u[posa+1][8,9,10]

	tmp3 = g_fuse(tmp2, 3)
	@grassmann tmp2[1,2,3,9,4,5,6,7;10] := tmp3[1,2,3,4,5,6,7,8] * v[posa+1][8,9,10]

	tmp3 = g_fuse(tmp2, 3)

	left = g_trace(tmp3, 2)
	# restore the accumulated left environment: single fermionic permute (via @grassmann)
	@grassmann lout[1; 2 3 4 5 6] := left[1,2,3,4,5,6]
	return lout
end

function update_pair_right(right::AbstractParityTensorMap{<:Number, 5, 1}, j::Int, x::Vector, y::Vector, z::Vector, u::Vector, v::Vector)
	posb = 2 * j

	@grassmann tmp1[4;5 6 7 1 2 8] := v[posb][1,2,3] * right[3,4,5,6,7,8]
	@grassmann tmp2[4;5 6 1 7 2 8 9] := u[posb][1,2,3] * tmp1[3,4,5,6,7,8,9]

	tmp3 = g_fuse(tmp2, 6)
	@grassmann tmp2[4;5 1 6 7 2 8 9] := z[posb][1,2,3] * tmp3[3,4,5,6,7,8,9]

	tmp3 = g_fuse(tmp2, 6)

	@grassmann tmp2[4;1 5 6 7 2 8 9] := y[posb][1,2,3] * tmp3[3,4,5,6,7,8,9]

	tmp3 = g_fuse(tmp2, 6)
	@grassmann tmp2[1;4 5 6 7 2 8 9] := x[posb][1,2,3] * tmp3[3,4,5,6,7,8,9]

	tmp3 = g_fuse(tmp2, 6)
	@grassmann tmp1[4; 1 2 5 6 7 8 9] := x[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9]

	@grassmann tmp2[6;1 4 2 5 7 8 9 10] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7,8,9,10]
	tmp3 = g_fuse(tmp2, 4)

	@grassmann tmp2[7;1 4 5 2 6 8 9 10] := z[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9,10]

	tmp3 = g_fuse(tmp2, 5)
	@grassmann tmp2[8;1 4 5 6 2 7 9 10] := u[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9,10]

	tmp3 = g_fuse(tmp2, 6)
	@grassmann tmp2[1;4 5 6 7 2 8 9 10] := v[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9,10]
	tmp3 = g_fuse(tmp2, 6)

	right = g_trace(tmp3, 6)
	# restore the accumulated right environment: single fermionic permute (via @grassmann)
	@grassmann rout[1 2 3 4 5; 6] := right[1,2,3,4,5,6]
	return rout
end
