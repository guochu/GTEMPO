function partialmpo_retardedinteract(row::Tuple{Int, Int}, cols::Vector{Tuple{Int, Int}}, coefs::Vector{<:Number}; trunc::TruncationScheme=DefaultMPOTruncation)
	# println("row=", row, " cols ", cols)
	@boundscheck begin
		(length(cols) == length(coefs) >= 1) || throw(ArgumentError("coefs and cols size mismatch"))
		issorted(cols) || throw(ArgumentError("cols should be sorted"))
		for item in cols
			x, y = item
			(x < y) || throw(ArgumentError("each item of cols should be sorted"))
		end
		for i in 1:length(cols)-1
			(cols[i][2] < cols[i+1][1]) || throw(ArgumentError("cols should be sorted"))
		end
	end
	# I2 = one(JW)

	ph = grassmannpspace()
	f = isomorphism(eltype(coefs), fuse(ph, ph), ph ⊗ ph)
	@tensor abar_a[1,7;9,8] := σ₊[1,2,3,4] * σ₋'[3,5,9,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor I4[5,6] := I2[1,2] * I2[3,4] * f[5,1,3] * conj(f[6,2,4])
    virtual = isomorphism(eltype(coefs), Rep[ℤ₂](0=>1), Rep[ℤ₂](0=>1))
    T = scalartype(virtual)
    @tensor m22I[1,3;2,4] := virtual[1,2] * I4[3,4] 
    # println(space_l(abar_a), " ", space_r(abar_a))
    ph2 = fuse(ph, ph)

    function leftmat(c) 
    	mat = Matrix{Any}(undef, 2, 2)
    	mat[1,1] = m22I
    	mat[1,2] = c * abar_a
   		mat[2,1] = 0
   		mat[2,2] = m22I
   		return SparseMPOTensor(mat, T, ph2)
    end
    function rightmat(c)
    	mat = Matrix{Any}(undef, 2, 2)
    	mat[1,1] = m22I
    	mat[1,2] = 0
    	mat[2,1] = c * abar_a
    	mat[2,2] = m22I
    	return SparseMPOTensor(mat, T, ph2)
    end
    function middlemat(c)
    	mat = Matrix{Any}(undef, 2, 2)
    	mat[1,1] = m22I + c * abar_a
    	# mat[1,1] = m22I 
    	mat[1,2] = abar_a
    	mat[2,1] = abar_a
    	mat[2,2] = 0 
    	return SparseMPOTensor(mat, T, ph2)
    end

	row_pos = findfirst(x -> row == x, cols)
	# isnothing(row_pos) && throw(ArgumentError("$(coefs) should contain $(row)"))

	positions2 = Int[]
	m = middlemat(0)
	if isnothing(row_pos)
		if row < cols[1]
			data = vcat(middlemat(0), [rightmat(coefs[i]) for i in 1:length(coefs)])
			push!(positions2, row[1])
			push!(positions2, row[2])
			for i in 1:length(coefs)
				push!(positions2, cols[i][1])
				push!(positions2, cols[i][2])
			end
		elseif row > cols[end]
			data = vcat([leftmat(coefs[i]) for i in 1:length(coefs)], middlemat(0))
			for i in 1:length(coefs)
				push!(positions2, cols[i][1])
				push!(positions2, cols[i][2])
			end	
			push!(positions2, row[1])
			push!(positions2, row[2])		
		else
			row_pos = findfirst(x -> row < x, cols)
			data = vcat([leftmat(coefs[i]) for i in 1:row_pos-1], middlemat(0), [rightmat(coefs[i]) for i in row_pos:length(coefs)])

			for i in 1:row_pos-1
				push!(positions2, cols[i][1])
				push!(positions2, cols[i][2])
			end
			push!(positions2, row[1])
			push!(positions2, row[2])
			for i in row_pos:length(coefs)
				push!(positions2, cols[i][1])
				push!(positions2, cols[i][2])
			end		
		end
	else
		if length(coefs) == 1
			data = [middlemat(coefs[1])]
		else
		    if row_pos == 1
		    	data = vcat(middlemat(coefs[1]), [rightmat(coefs[i]) for i in 2:length(coefs)])
		    elseif row_pos == length(cols)
		    	data = vcat([leftmat(coefs[i]) for i in 1:length(coefs)-1], middlemat(coefs[end]))
		    else
		    	data = vcat([leftmat(coefs[i]) for i in 1:row_pos-1], middlemat(coefs[row_pos]), [rightmat(coefs[i]) for i in row_pos+1:length(coefs)])
		    end
		end
	  
	    for item in cols
	    	push!(positions2, item[1])
	    	push!(positions2, item[2])
	    end

	end

	mpo = MPO(MPOHamiltonian(data))
	data2 = similar(mpo.data, 2*length(data))
	for i in 1:length(mpo)
		data2[2*i-1], data2[2*i] = split_mpotensor(mpo[i], trunc)
	end
	return PartialMPO(data2, positions2)
end