

"""
	updateright(hold::AbstractArray{T, 2}, mpsAj::MPSTensor, mpsBj::MPSTensor) where {T<:Number}
update storage from right to left for overlap of mps
"""
function updateright(hold::AbstractMatrix, mpsAj::DenseMPSTensor, mpsBj::DenseMPSTensor) 
	@tensor m2[-1 -2;-3] := conj(mpsAj[-1, -2, 1]) * hold[1, -3]
	@tensor hnew[-1;-2] := m2[-1,1,2] * mpsBj[-2,1,2]
	return hnew
end

"""
	updateleft(hold::AbstractArray{T, 2}, mpsAj::MPSTensor, mpsBj::MPSTensor) where {T<:Number}
update storage from left to right for overlap of mps
"""
function updateleft(hold::AbstractMatrix, mpsAj::DenseMPSTensor, mpsBj::DenseMPSTensor) 
	@tensor m2[-1 -2; -3] := conj(mpsAj[1, -2, -3]) * hold[1, -1]
	@tensor hnew[-1; -2] := m2[1,2,-1] * mpsBj[1,2,-2]
	return hnew
end


function n_fuse(m::AbstractArray{<:Number, N}, i::Int) where {N}
	@assert (1 <= i < N)
	@assert size(m, i) == size(m, i+1) == 2
	# @assert (i < M) || (M <= i < N)
	local tmp
	s_front = ntuple(x->size(m, x), i-1)
	s_tail = ntuple(x->size(m, i+1+x), N-i-1)
	s1 = prod(s_front)
	s2 = prod(s_tail)

	m4 = reshape(m, s1, 2, 2, s2)
	m3 = zeros(s1, 2, s2)
	m3[:,1,:] = m4[:, 1, 1, :]
	m3[:,2,:] = m4[:, 1, 2, :] + m4[:, 2, 1, :] + m4[:, 2, 2, :] 

	new_size = (s_front..., 2, s_tail...)
	return reshape(m3, new_size)
end