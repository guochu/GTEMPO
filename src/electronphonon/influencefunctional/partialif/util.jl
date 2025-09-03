

# function partialdensempo(row::Int, cols::Vector{Int}, coefs::Vector{<:Number})
# 	# println("row=", row, " cols ", cols)
# 	@assert length(cols) == length(coefs) > 1
# 	# @assert issorted(cols)
# 	p = sortperm(cols)
# 	cols = cols[p]
# 	coefs = coefs[p]
# 	T = scalartype(coefs)
# 	I2 = _eye(T, 2)
# 	_zero = zeros(T, 2, 2)
# 	nop = Matrix{T}([0 0; 1 1])

# 	function leftboundaryblock(m) 
# 		tmp = Matrix{Matrix{T}}(undef, 1, 2)
# 		tmp[1,1] = I2
# 		tmp[1,2] = m
# 		return tmp
# 	end
# 	function rightboundaryblock(m)
# 		tmp = Matrix{Matrix{T}}(undef, 2, 1)
# 		tmp[1,1] = I2
# 		tmp[2,1] = m
# 		return tmp
# 	end
# 	function leftblock(m)
# 		tmp = Matrix{Matrix{T}}(undef, 2, 2)
#         tmp[1,1] = I2
#         tmp[2,2] = I2
#         tmp[1,2] = _zero
#         tmp[2,1] = m
#         return tmp
# 	end
# 	function rightblock(m)
#         tmp = Matrix{Matrix{T}}(undef, 2, 2)
#         tmp[1,1] = I2
#         tmp[2,2] = I2
#         tmp[2,1] = _zero
#         tmp[1,2] = m
#         return tmp	
# 	end

# 	if row < cols[1]
# 	    data = [leftboundaryblock(nop)]
# 	    for i in 1:length(cols)-1
# 	        push!(data, leftblock(nop * coefs[i]))
# 	    end
# 	    push!(data, rightboundaryblock(nop * coefs[end]))
# 	    positions = [row; cols]
# 	elseif row == cols[1]
# 	    tmp = Matrix{Matrix{T}}(undef, 1, 2)
# 	    tmp[1, 1] = I2 + coefs[1] * nop
# 	    tmp[1, 2] = nop	
# 	    data = [tmp]
# 	    for i in 2:length(cols)-1
# 	        push!(data, leftblock(nop * coefs[i]))
# 	    end
# 	    push!(data, rightboundaryblock(nop * coefs[end]))
# 	    positions = cols
# 	elseif row > cols[end]
# 	    data = [leftboundaryblock(coefs[1] * nop)]
# 	    for i in 2:length(cols)
# 	        push!(data, rightblock(coefs[i] * nop))
# 	    end	
# 	    push!(data, rightboundaryblock(nop))
# 	    positions = [cols; row]
# 	elseif row == cols[end]
# 	    data = [leftboundaryblock(coefs[1] * nop)]
# 	    for i in 2:length(cols)-1
# 	        push!(data, rightblock(coefs[i] * nop))
# 	    end	
# 	    tmp = Matrix{Matrix{T}}(undef, 2, 1)
# 	    tmp[1, 1] = I2 + coefs[end] * nop
# 	    tmp[2, 1] = nop	
# 	    push!(data, tmp)
# 	    positions = cols
# 	else
# 		row_pos = findfirst(x -> row == x, cols)
# 		if isnothing(row_pos)
# 			data = [leftboundaryblock(coefs[1] * nop)]
# 			for i in 2:row_pos-1
# 		        push!(data, leftblock(coefs[i] * nop))			
# 			end			
# 			tmp = Matrix{Matrix{T}}(undef, 2, 2)
# 			tmp[1,1] = I2
# 			tmp[2,2] = _zero
# 			tmp[1,2] = nop
# 			tmp[2,1] = nop
# 			push!(data, tmp)

# 		    for i in row_pos:length(cols)-1
# 		        push!(data, rightblock(coefs[i]*nop))
# 		    end
# 		    push!(data, rightboundaryblock(coefs[end]*nop))
# 		    positions = insert!(copy(cols), row_pos, row)
# 		else
# 			data = [leftboundaryblock(coefs[1] * nop)]
# 			for i in 2:row_pos-1
# 		        push!(data, leftblock(coefs[i] * nop))			
# 			end			
# 			tmp = Matrix{Matrix{T}}(undef, 2, 2)
# 			tmp[1,1] = I2 + coefs[row_pos] * nop
# 			tmp[2,2] = _zero
# 			tmp[1,2] = nop
# 			tmp[2,1] = nop
# 			push!(data, tmp)

# 		    for i in row_pos+1:length(cols)-1
# 		        push!(data, rightblock(coefs[i]*nop))
# 		    end
# 		    push!(data, rightboundaryblock(coefs[end]*nop))
# 		    positions = cols
# 		end
# 	end
# 	mpodata = _convert_to_mpotensor.(data)
# 	return PartialDenseMPO(mpodata, positions)
# end


# the algorithm in StrathearnLovett2018
function partialif_densemps(lat_size::Int, row::Int, cols::Vector{Int}, coefs::Vector{<:Number})
	# println("row=", row, " cols ", cols)
	@assert length(cols) == length(coefs) > 1
	# @assert (row in cols)
	# @assert issorted(cols)
	p = sortperm(cols)
	cols = cols[p]
	coefs = exp.(coefs[p])
	T = scalartype(coefs)

	function onebody(m) 
		tmp = ones(T,2)
		tmp[1] = m
		return tmp
	end
	function twobody(m)
		tmp = ones(T, 2,2)
		tmp[1,1] = m
		return tmp
	end
	L = length(cols)
	mpsdata = Vector{Array{T, 3}}(undef, L)
	pos = findfirst(x->row==x, cols)
	isnothing(pos) && throw(ArgumentError("$(row) is not a member of $cols"))
	if pos == 1
		m = onebody(coefs[1])
		tmp = zeros(T, 2, 2)
		for i in 1:2
			tmp[i, i] = m[i]
		end
		mpsdata[1] = reshape(tmp, 1,2,2)
		for j in 2:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, 2,2,2)
			for i in 1:2
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp
		end
		m = twobody(coefs[L])
		mpsdata[L] = reshape(m,2,2,1)
	elseif pos == L
		m = twobody(coefs[1])
		mpsdata[1] = reshape(m, 1,2,2)
		for j in 2:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, 2,2,2)
			for i in 1:2
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp			
		end
		m = onebody(coefs[L])
		tmp = zeros(T, 2, 2)
		for i in 1:2
			tmp[i, i] = m[i]
		end
		mpsdata[L] = reshape(tmp,2,2,1)
	else
		m = twobody(coefs[1])
		mpsdata[1] = reshape(m, 1,2,2)
		for j in 2:pos-1
			m = twobody(coefs[j])
			tmp = zeros(T, 2,2,2)
			for i in 1:2
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp
		end
		m = onebody(coefs[pos])
		tmp = zeros(T,2,2,2)
		for i in 1:2
			tmp[i,i,i] = m[i]
		end
		mpsdata[pos] = tmp
		for j in pos+1:L-1
			m = twobody(coefs[j])
			tmp = zeros(T, 2,2,2)
			for i in 1:2
				tmp[i,:,i] = m[i,:]
			end
			mpsdata[j] = tmp			
		end
		m = twobody(coefs[L])
		mpsdata[L] = reshape(m,2,2,1)
	end
	# return mpsdata
	return _fit_to_full(lat_size, cols, mpsdata)
end

function _fit_to_full(L::Int, pos, mpsdata)
	r = similar(mpsdata, L)
	for j in 1:pos[1]-1
		r[j] = ones(1,2,1)
	end
	for j in pos[end]+1:L
		r[j] = ones(1,2,1)
	end
	tmp = zeros(2,2,2)
	for i in 1:2
		tmp[:,i,:] = _eye(2)
	end
	for j in pos[1]:pos[end]
		posj = findfirst(x->x==j, pos)
		r[j] = isnothing(posj) ? tmp : mpsdata[posj]
	end
	return FockMPS(r)
end

# function _convert_to_mpotensor(data::Matrix{Matrix{T}}) where {T<:Number}
# 	s1, s2 = size(data)
# 	d1, d2 = size(data[1])
# 	r = zeros(T, s1, d1, s2, d2)
# 	for i in 1:s1, j in 1:s2
# 		r[i, :, j, :] = data[i, j]
# 	end
# 	return r
# end
