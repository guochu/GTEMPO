

# function partialdensempo(row::Int, cols::Vector{Int}, coefs::Vector{<:Number})
# 	# println("row=", row, " cols ", cols)
# 	@assert length(cols) == length(coefs) >= 1
# 	# @assert issorted(cols)
# 	p = sortperm(cols)
# 	cols = cols[p]
# 	coefs = coefs[p]
# 	T = scalartype(coefs)
# 	I2 = _eye(T, 2)
# 	_zero = zeros(T, 2, 2)
# 	nop = Matrix{scalartype(x)}([0 0; 1 1])

# 	if row < cols[1]
# 	    tmp = Matrix{Matrix{T}}(undef, 1, 2)
# 	    tmp[1, 1] = I2
# 	    tmp[1, 2] = nop
# 	    data = [tmp]
# 	    for i in 1:length(cols)-1
# 	        tmp = Matrix{Matrix{T}}(undef, 2, 2)
# 	        tmp[1,1] = I2
# 	        tmp[2,2] = I2
# 	        tmp[1,2] = _zero
# 	        tmp[2,1] = nop * coefs[i]
# 	        push!(data, tmp)
# 	    end
# 	    tmp = Matrix{Matrix{T}}(undef, 2, 1)
# 	    tmp[1,1] = I2
# 	    tmp[2,1] = nop * coefs[end]
# 	    push!(data, tmp)
# 	    positions = [row; cols]
# 	elseif row > cols[end]
# 	    tmp = Matrix{Matrix{T}}(undef, 1, 2)
# 	    tmp[1, 1] = I2
# 	    tmp[1, 2] = coefs[1] * nop
# 	    data = [tmp]
# 	    for i in 2:length(cols)
# 	        tmp = Matrix{Matrix{T}}(undef, 2, 2)
# 	        tmp[1,1] = I2
# 	        tmp[2,2] = I2
# 	        tmp[2,1] = _zero
# 	        tmp[1,2] = coefs[i] * nop
# 	        push!(data, tmp)
# 	    end	
# 	    tmp = Matrix{Matrix{T}}(undef, 2, 1)
# 	    tmp[1,1] = I2
# 	    tmp[2,1] = nop
# 	    push!(data, tmp)
# 	    positions = [cols; row]
# 	else
# 	    tmp = Matrix{Matrix{T}}(undef, 1, 2)
# 	    tmp[1, 1] = I2
# 	    tmp[1, 2] = coefs[1] * nop
# 	    data = [tmp]

# 		row_pos = findfirst(x -> row < x, cols)
# 		# println("row pos is ", row_pos)
# 		for i in 2:row_pos-1
# 	        tmp = Matrix{Any}(undef, 2, 2)
# 	        tmp[1,1] = I2
# 	        tmp[2,2] = m22I
# 	        tmp[2,1] = 0.
# 	        tmp[1,2] = (-coefs[i]) * σ₊
# 	        push!(data, SparseMPOTensor(tmp, T, pspace))			
# 		end
# 		tmp = Matrix{Any}(undef, 2, 2)
# 		tmp[1,1] = I2
# 		tmp[2,2] = 0.
# 		tmp[1,2] = σ₊
# 		tmp[2,1] = adjoint(σ₋) * one(eltype(coefs))
# 		push!(data, SparseMPOTensor(tmp, T, pspace))

# 	    for i in row_pos:length(cols)-1
# 	        tmp = Matrix{Any}(undef, 2, 2)
# 	        tmp[1,1] = I2
# 	        tmp[2,2] = m22I
# 	        tmp[1,2] = 0.
# 	        tmp[2,1] = adjoint(σ₋) * coefs[i]
# 	        push!(data, SparseMPOTensor(tmp, T, pspace))
# 	    end
# 	    tmp = Matrix{Any}(undef, 2, 1)
# 	    tmp[1,1] = I2
# 	    tmp[2,1] = adjoint(σ₋) * coefs[end]
# 	    mpof = SparseMPOTensor(tmp, T, pspace)
# 	    push!(data, mpof)
# 	    positions = insert!(copy(cols), row_pos, row)
# 	end
# 	mpo = MPO(MPOHamiltonian(data))
# 	return PartialMPO(mpo.data, positions)
# end



# function split_mpotensor(mpoj::MPOTensor, trunc)
# 	ph = grassmannpspace()
# 	f = isomorphism(scalartype(mpoj), fuse(ph, ph), ph ⊗ ph)
# 	@tensor mpoj6[1,5,7;6,3,8] := mpoj[1,2,3,4] * conj(f[2,5,6]) * f[4,7,8]
# 	u, s, v = tsvd!(mpoj6, trunc=trunc)
# 	ss = sqrt(s)
# 	return permute(u * ss, (1,2), (4,3)), permute(ss * v, (1,2), (3,4))
# end


