

# function influenceoperator(lattice::ImagGrassmannLattice1Order{<:A1A1B1B1}, corr::ImagCorrelationFunction; band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
# 	mpoj = ti_mpotensor(corr, algexpan)
# 	h = MPOHamiltonian([mpoj, mpoj, mpoj])
# 	mpo = MPO(h)
# 	# println("here********", space_l(mpo), " ", space_r(mpo))
# 	u_left, v_left = _split_op(mpo[1], DefaultMPOTruncation)
# 	u_middle, v_middle = _split_op(mpo[2], DefaultMPOTruncation)
# 	u_right, v_right = _split_op(mpo[3], DefaultMPOTruncation)
# 	L = length(lattice)
# 	data2 = Vector{typeof(u_middle)}(undef, L)
# 	leftspace = oneunit(grassmannpspace())
# 	_JW = JW
# 	# for band in 1:lattice.bands
# 	# 	for c in (true, false)
# 	# 		pos = index(lattice, 0, band=band, conj=c)
# 	# 		vd2 = id(leftspace)
# 	# 		data2[pos] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 	# 	end
# 	# end
# 	for pos in 1:index(lattice, lattice.k, conj=true, band=lattice.bands)
# 		vd2 = id(leftspace)
# 		data2[pos] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 	end
# 	j = lattice.N
# 	posa, posb = band_boundary(lattice, j)
# 	pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
# 	@assert (posa <= pos1 < pos2 <= posb) 
# 	for i in posa:posb
# 		if i < pos1
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 		elseif i == pos1
# 			data2[i] = u_left
# 		elseif i < pos2
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
# 		elseif i == pos2
# 			data2[i] = v_left
# 		else
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
# 		end
# 		leftspace = space_r(data2[i])'
# 	end
# 	for j in lattice.N-1:-1:2
# 		posa, posb = band_boundary(lattice, j)
# 		pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
# 		@assert (posa <= pos1 < pos2 <= posb) 
# 		for i in posa:posb
# 			if i < pos1
# 				vd2 = id(leftspace)
# 				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 			elseif i == pos1
# 				data2[i] = u_middle
# 			elseif i < pos2
# 				vd2 = id(leftspace)
# 				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
# 			elseif i == pos2
# 				data2[i] = v_middle
# 			else
# 				vd2 = id(leftspace)
# 				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
# 			end
# 			leftspace = space_r(data2[i])'
# 		end
# 	end
# 	j = 1
# 	posa, posb = band_boundary(lattice, j)
# 	pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
# 	@assert (posa <= pos1 < pos2 <= posb) 
# 	for i in posa:posb
# 		if i < pos1
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 		elseif i == pos1
# 			data2[i] = u_right
# 		elseif i < pos2
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
# 		elseif i == pos2
# 			data2[i] = v_right
# 		else
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
# 		end
# 		leftspace = space_r(data2[i])'
# 	end
# 	return MPO(data2)
# end

# function differentialinfluenceoperator(lattice::ImagGrassmannLattice1Order{<:A1A1B1B1}, corr::ImagCorrelationFunction, dt::Real, alg::TimeEvoMPOAlgorithm=WII(); 
# 										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
# 	mpoj = ti_mpotensor(corr, algexpan)
# 	h = MPOHamiltonian([mpoj, mpoj, mpoj])
# 	h2 = timeevompo(h, dt, alg)
# 	mpo = MPO(h2)
# 	# println("here********", space_l(mpo), " ", space_r(mpo))
# 	u_left, v_left = _split_op(mpo[1], DefaultMPOTruncation)
# 	u_middle, v_middle = _split_op(mpo[2], DefaultMPOTruncation)
# 	u_right, v_right = _split_op(mpo[3], DefaultMPOTruncation)
# 	L = length(lattice)
# 	data2 = Vector{typeof(u_middle)}(undef, L)
# 	leftspace = oneunit(grassmannpspace())
# 	_JW = I2
# 	for pos in 1:index(lattice, lattice.k, conj=true, band=lattice.bands)
# 		vd2 = id(leftspace)
# 		data2[pos] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 	end
# 	j = lattice.N
# 	posa, posb = band_boundary(lattice, j)
# 	pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
# 	@assert (posa <= pos1 < pos2 <= posb) 
# 	for i in posa:posb
# 		if i < pos1
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 		elseif i == pos1
# 			data2[i] = u_left
# 		elseif i < pos2
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
# 		elseif i == pos2
# 			data2[i] = v_left
# 		else
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
# 		end
# 		leftspace = space_r(data2[i])'
# 	end
# 	for j in lattice.N-1:-1:2
# 		posa, posb = band_boundary(lattice, j)
# 		pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
# 		@assert (posa <= pos1 < pos2 <= posb) 
# 		for i in posa:posb
# 			if i < pos1
# 				vd2 = id(leftspace)
# 				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 			elseif i == pos1
# 				data2[i] = u_middle
# 			elseif i < pos2
# 				vd2 = id(leftspace)
# 				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
# 			elseif i == pos2
# 				data2[i] = v_middle
# 			else
# 				vd2 = id(leftspace)
# 				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
# 			end
# 			leftspace = space_r(data2[i])'
# 		end
# 	end
# 	j = 1
# 	posa, posb = band_boundary(lattice, j)
# 	pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
# 	@assert (posa <= pos1 < pos2 <= posb) 
# 	for i in posa:posb
# 		if i < pos1
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
# 		elseif i == pos1
# 			data2[i] = u_right
# 		elseif i < pos2
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
# 		elseif i == pos2
# 			data2[i] = v_right
# 		else
# 			vd2 = id(leftspace)
# 			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
# 		end
# 		leftspace = space_r(data2[i])'
# 	end
# 	return MPO(data2)
# end



# build the translationally invariant MPOTensor of IF
function ti_mpotensor(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm)
	ph = grassmannpspace()
	f = isomorphism(fuse(ph, ph), ph ⊗ ph)
	_JW = JW
	@tensor a[1,7;3,8] := σ₊[1,2,3,4] * _JW[5,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor abar[1,7;3,8] := _JW[2,4] * σ₋'[1,5,3,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor JW4[5;6] := _JW[1,2] * _JW[3,4] * f[5,1,3] * conj(f[6,2,4])
	m1 = GenericDecayTerm(a, abar, -corr.ηₖⱼ[2:end], middle = JW4)

	@tensor a[1,7;3,8] := σ₋'[1,2,3,4] * I2[5,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor abar[1,7;3,8] := I2[2,4] * σ₊[1,5,3,6] * f[7,2,5] * conj(f[8,4,6])
	m2 = GenericDecayTerm(abar, a, corr.ηⱼₖ[2:end], middle = JW4)

	m1s = exponential_expansion(m1, alg=alg)
	m2s = exponential_expansion(m2, alg=alg)
	@tensor abar_a[7;8] := σ₊[1,2,3,4] * σ₋'[3,5,1,6] * f[7,2,5] * conj(f[8,4,6])

	coef = corr.ηⱼₖ[1] + corr.ηₖⱼ[1]
	return SchurMPOTensor(-coef * abar_a, vcat(m1s, m2s))
end

function _split_op(mpoj::MPOTensor, trunc)
	ph = grassmannpspace()
	f = isomorphism(fuse(ph, ph), ph ⊗ ph)
	@tensor mpoj6[1,5,7;6,3,8] := mpoj[1,2,3,4] * conj(f[2,5,6]) * f[4,7,8]
	u, s, v = tsvd!(mpoj6, trunc=trunc)
	ss = sqrt(s)
	return permute(u * ss, (1,2), (4,3)), permute(ss * v, (1,2), (3,4))
end

