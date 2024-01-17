# 2 gmps integration

function update_pair_left(left::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS; trunc=DefaultIntegrationTruncation)
	f = (scaling(x) * scaling(y))^2
	posa = 2*j-1
	@tensor twositemps1[4,3,5;6] := left[1,2]*x[posa][1,3,4]*y[posa][2,5,6]
	mpsj1 = _fuse_physical(swap12!(twositemps1))
	m1, m2 = _swap_gate(mpsj1, y[posa+1], trunc=trunc)
	@tensor twositemps2[3,2,4;5] := x[posa+1][1,2,3] * m1[1,4,5]
	mpsj2 = _fuse_physical(swap12!(twositemps2))
	@tensor twositemps3[1,2,4;5] := mpsj2[1,2,3] * m2[3,4,5]
	left = _trace_physical(twositemps3, nt=false)
	left = rmul!(left, f)
	return left	
end

function update_pair_right(right::AbstractTensorMap, j::Int, x::GrassmannMPS, y::GrassmannMPS; trunc=DefaultIntegrationTruncation)
	f = (scaling(x) * scaling(y))^2
	posb = 2 * j
	@tensor twositemps1[1; 2 6 5] := y[posb][1,2,3] * right[3,4] *x[posb][5,6,4]
	mpsj1 = _fuse_physical(swap34!(twositemps1))
	m1, m2 = _swap_gate(y[posb-1], mpsj1, trunc=trunc)
	@tensor twositemps2[1; 2 5 4] := m2[1,2,3] * x[posb-1][4,5,3] 
	mpsj2 = _fuse_physical(swap34!(twositemps2))
	@tensor twositemps3[1,2,4;5] := m1[1,2,3] * mpsj2[3,4,5] 
	right = _trace_physical(twositemps3, nt=false)
	right = rmul!(right, f )
	return right
end

