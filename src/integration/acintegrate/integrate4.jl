
function update_pair_left(left::AbstractParityTensorMap{<:Number, 1, 4}, j::Int, x::Vector, y::Vector, z::Vector, u::Vector)
	posa = 2*j-1

	@grassmann tmp1[7,1,2,5,6; 3] := left[7,1,2,3,4] * u[posa][4,5,6]
	@grassmann tmp2[8,1,3,6,4,7;2] := tmp1[8,1,2,3,4,5] * z[posa][5,6,7]

	tmp3 = g_fuse(tmp2, 3)

	@grassmann tmp2[8,2,6,3,4,7; 1] := tmp3[8,1,2,3,4,5] * y[posa][5,6,7]

	tmp3 = g_fuse(tmp2, 2)

	@grassmann tmp2[8,1,6,2,3,4;7] := tmp3[8,1,2,3,4,5] * x[posa][5,6,7]

	tmp3 = g_fuse(tmp2, 2)

	# \bar{a}
	@grassmann tmp1[8,1,2,3,6,7;4] := tmp3[8,1,2,3,4,5] * x[posa+1][5,6,7]
	@grassmann tmp2[9,1,2,4,7,5,8;3] := tmp1[9,1,2,3,4,5,6] * y[posa+1][6,7,8]

	# fuse indices
	tmp3 = g_fuse(tmp2, 4)

	@grassmann tmp2[9,1,3,7,4,5,8;2] := tmp3[9,1,2,3,4,5,6] * z[posa+1][6,7,8]

	# fuse indices
	tmp3 = g_fuse(tmp2, 3)

	@grassmann tmp2[9,1,2,7,3,4;5,8] := tmp3[9,1,2,3,4,5,6] * u[posa+1][6,7,8]
	# fuse indices
	tmp3 = g_fuse(tmp2, 3)

	# trace physices
	left = g_trace(tmp3, 2)
	# restore the accumulated left environment: single fermionic permute (via @grassmann)
	@grassmann lout[1; 2 3 4 5] := left[1,2,3,4,5]
	return lout
end


function update_pair_right(right::AbstractParityTensorMap{<:Number, 4, 1}, j::Int, x::Vector, y::Vector, z::Vector, u::Vector)
	posb = 2 * j
	@grassmann tmp1[4 ;5 6 1 2 7] := u[posb][1,2,3] * right[3,4,5,6,7]
	@grassmann tmp2[4;5 1 6 2 7 8] := z[posb][1,2,3] * tmp1[3,4,5,6,7,8]

	# fuse indices
	tmp3 = g_fuse(tmp2, 5)
	# println(space(y[posb], 3), " ", space(tmp2,1), " ", space(tmp3,1), " ", space(right, 3))

	@grassmann tmp2[4;1 5 6 2 7 8] := y[posb][1,2,3] * tmp3[3,4,5,6,7,8]

	# fuse indices
	tmp3 = g_fuse(tmp2, 5)

	@grassmann tmp2[1 ;4 5 6 2 7 8] := x[posb][1,2,3] * tmp3[3,4,5,6,7,8]

	# fuse indices
	tmp3 = g_fuse(tmp2, 5)

	@grassmann tmp1[4 ;5 6 1 2 7 8] := x[posb-1][1,2,3] * tmp3[3,4,5,6,7,8]
	@grassmann tmp2[4;5 1 6 2 7 8 9] := y[posb-1][1,2,3] * tmp1[3,4,5,6,7,8,9]

	tmp3 = g_fuse(tmp2, 5)

	@grassmann tmp2[4;1 5 6 2 7 8 9] := z[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9]
	tmp3 = g_fuse(tmp2, 5)

	@grassmann tmp2[1 ;4 5 6 2 7 8 9] := u[posb-1][1,2,3] * tmp3[3,4,5,6,7,8,9]
	tmp3 = g_fuse(tmp2, 5)

	right = g_trace(tmp3, 5)
	# restore the accumulated right environment: single fermionic permute (via @grassmann)
	@grassmann rout[1 2 3 4; 5] := right[1,2,3,4,5]
	return rout

end
