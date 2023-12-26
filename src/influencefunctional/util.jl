# move left and right most physical indices after the i-th position
function pif_util(data::Vector{M}, i::Int; trunc) where {M}
	L = length(data) 
	ndata = Vector{M}(undef, L+1)
	util = TensorMap(ds->ones(Float64, ds), oneunit(_ph) ← one(_ph))
	@tensor left[1,3;4,5,2] := util[1] * data[1][2,3,4,5]
	for j in 1:i-1
		ndata[j], left = move_to_right(left, data[j+1], trunc=trunc)
	end
	if i == L
		u, s, v = DMRG.stable_tsvd(left, (1,2,4), (3,5), trunc=trunc)
		ndata[i] = permute(u, (1,2), (4,3))
		ndata[i+1] = @tensor tmp[1,4;5,3] := s[1,2] * v[2,3,4] * conj(util[5])
	else
		@tensor right[3,1,2;5,4] := data[L][1,2,3,4] * conj(util[5])
		for j in (L+1):-1:(i+3)
			right, ndata[j] = move_to_left(data[j-2], right, trunc=trunc)
		end
		if i == 0
			u, s, v = DMRG.stable_tsvd(right, (1,2), (3,4,5), trunc=trunc)
			u = u * s
			ndata[2] = permute(v, (1,2), (3,4))
			ndata[1] = @tensor tmp[1,3;4,2] := util[1] * u[2,3,4]
		else
			ndata[i], ndata[i+1], ndata[i+2] = central(left, right, trunc=trunc)	
		end			
	end
	return ndata
end

# t1 is rank-4, t2 is rank-5
function move_to_left(t1, t2; trunc)
	@tensor tmp[5,1,2,4; 6,7,8] := t1[1,2,3,4] * t2[5,3,6,7,8]
	u, s, v, err = DMRG.stable_tsvd!(tmp, trunc=trunc)
	return permute(u*s, (1,2,3), (5,4)), permute(v, (1,2), (3,4))
end

# t1 is rank-5, t2 is rank-4
function move_to_right(t1, t2; trunc)
	@tensor tmp[1,2,4; 6,7,8,5] := t1[1,2,3,4,5] * t2[3,6,7,8]
	u, s, v, err = DMRG.stable_tsvd!(tmp, trunc=trunc)
	return permute(u, (1,2), (4,3)), permute(s * v, (1,2), (3,4,5))
end

# t1 is rank-5, t2 is rank-5
function central(t1, t2; trunc)
	u1, s1, v1, err = DMRG.stable_tsvd(t1, (1,2,4), (3,5), trunc=trunc)
	v1 = s1 * v1
	u2, s2, v2, err = DMRG.stable_tsvd(t2, (1,2), (3,4,5), trunc=trunc)
	u2 = u2 * s2
	@tensor center[1,3;5,4] := v1[1,2,3] * u2[4,2,5]
	return permute(u1, (1,2), (4,3)), center, permute(v2, (1,2), (3,4))
end


function partialmpo(row::Int, cols::Vector{Int}, coefs::Vector{<:Number})
	# println("row=", row, " cols ", cols)
	@assert length(cols) == length(coefs) >= 1
	# @assert issorted(cols)
	p = sortperm(cols)
	cols = cols[p]
	coefs = coefs[p]
	I2 = one(JW)

    virtual = isomorphism(Matrix{eltype(coefs)}, Rep[ℤ₂](1=>1), Rep[ℤ₂](1=>1))
    pspace = grassmannpspace()
    T = scalartype(virtual)
    @tensor m22I[1,3;2,4] := virtual[1,2] * I2[3,4] 
	if row < cols[1]
	    tmp = Matrix{Any}(undef, 1, 2)
	    tmp[1, 1] = I2
	    tmp[1, 2] = one(eltype(coefs)) * σ₊
	    mpoj = SparseMPOTensor(tmp, T, pspace)
	    data = [mpoj]
	    for i in 1:length(cols)-1
	        tmp = Matrix{Any}(undef, 2, 2)
	        tmp[1,1] = I2
	        tmp[2,2] = m22I
	        tmp[1,2] = 0.
	        tmp[2,1] = adjoint(σ₋) * coefs[i]
	        push!(data, SparseMPOTensor(tmp, T, pspace))
	    end
	    tmp = Matrix{Any}(undef, 2, 1)
	    tmp[1,1] = I2
	    tmp[2,1] = adjoint(σ₋) * coefs[end]
	    mpof = SparseMPOTensor(tmp, T, pspace)
	    push!(data, mpof)
	    positions = [row; cols]
	elseif row > cols[end]
	    tmp = Matrix{Any}(undef, 1, 2)
	    tmp[1, 1] = I2
	    tmp[1, 2] = (-coefs[1]) * σ₊
	    mpoj = SparseMPOTensor(tmp, T, pspace)
	    data = [mpoj]
	    for i in 2:length(cols)
	        tmp = Matrix{Any}(undef, 2, 2)
	        tmp[1,1] = I2
	        tmp[2,2] = m22I
	        tmp[2,1] = 0.
	        tmp[1,2] = (-coefs[i]) * σ₊
	        push!(data, SparseMPOTensor(tmp, T, pspace))
	    end	
	    tmp = Matrix{Any}(undef, 2, 1)
	    tmp[1,1] = I2
	    tmp[2,1] = adjoint(σ₋) * one(eltype(coefs))
	    mpof = SparseMPOTensor(tmp, T, pspace)
	    push!(data, mpof)
	    positions = [cols; row]
	else
	    tmp = Matrix{Any}(undef, 1, 2)
	    tmp[1, 1] = I2
	    tmp[1, 2] = (-coefs[1]) * σ₊
	    mpoj = SparseMPOTensor(tmp, T, pspace)
	    data = [mpoj]

		row_pos = findfirst(x -> row < x, cols)
		# println("row pos is ", row_pos)
		for i in 2:row_pos-1
	        tmp = Matrix{Any}(undef, 2, 2)
	        tmp[1,1] = I2
	        tmp[2,2] = m22I
	        tmp[2,1] = 0.
	        tmp[1,2] = (-coefs[i]) * σ₊
	        push!(data, SparseMPOTensor(tmp, T, pspace))			
		end
		tmp = Matrix{Any}(undef, 2, 2)
		tmp[1,1] = I2
		tmp[2,2] = 0.
		tmp[1,2] = σ₊
		tmp[2,1] = adjoint(σ₋) * one(eltype(coefs))
		push!(data, SparseMPOTensor(tmp, T, pspace))

	    for i in row_pos:length(cols)-1
	        tmp = Matrix{Any}(undef, 2, 2)
	        tmp[1,1] = I2
	        tmp[2,2] = m22I
	        tmp[1,2] = 0.
	        tmp[2,1] = adjoint(σ₋) * coefs[i]
	        push!(data, SparseMPOTensor(tmp, T, pspace))
	    end
	    tmp = Matrix{Any}(undef, 2, 1)
	    tmp[1,1] = I2
	    tmp[2,1] = adjoint(σ₋) * coefs[end]
	    mpof = SparseMPOTensor(tmp, T, pspace)
	    push!(data, mpof)
	    positions = insert!(copy(cols), row_pos, row)
	end
	mpo = MPO(MPOHamiltonian(data))
	return QTerm(positions, mpo.data)
end
