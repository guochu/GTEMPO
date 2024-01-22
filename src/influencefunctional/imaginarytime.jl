function partialinfluencefunctional(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1)
	row = index(lattice, i, band=band, conj=true)
	col_pos = [index(lattice, j, band=band, conj=false) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end

"""
	influenceoperator(lattice::ImagGrassmannLattice1Order{<:A1A1B1B1}, corr2::ImagCorrelationFunction; band, algexpan)
"""
function influenceoperator(lattice::ImagGrassmannLattice1Order{<:A1A1B1B1}, corr2::ImagCorrelationFunction; band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr.data
	mpoj = ti_mpotensor(corr, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	mpo = MPO(h)
	_JW = JW
	return _fit_to_lattice(lattice, mpo, _JW, band) 
end

"""
	differentialinfluencefunctional(lattice::ImagGrassmannLattice1Order{<:A1A1B1B1}, corr2::ImagCorrelationFunction, dt, alg; band, algexpan)
"""
function differentialinfluencefunctional(lattice::ImagGrassmannLattice1Order{<:A1A1B1B1}, corr2::ImagCorrelationFunction, dt::Real, alg::FirstOrderStepper; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr.data
	mpoj = ti_mpotensor(corr, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	h2 = timeevompo(h, dt, alg)
	mpo = MPO(h2)
	_JW = I2
	return _fit_to_lattice(lattice, mpo, _JW, band) 
end

function _fit_to_lattice(lattice::ImagGrassmannLattice, mpo::MPO, _JW::MPSBondTensor, band::Int)
	@assert length(mpo) == 3
	(LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("only works lattice with TimeLocalLayout"))

	u_left, v_left = _split_op(mpo[1], DefaultMPOTruncation)
	u_middle, v_middle = _split_op(mpo[2], DefaultMPOTruncation)
	u_right, v_right = _split_op(mpo[3], DefaultMPOTruncation)
	L = length(lattice)
	data2 = Vector{typeof(u_middle)}(undef, L)
	leftspace = oneunit(grassmannpspace())

	posa, posb = band_boundary(lattice, 0)
	for pos in posa:posb
		vd2 = id(leftspace)
		data2[pos] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
	end
	j = lattice.k
	posa, posb = band_boundary(lattice, j)
	pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
	@assert (posa <= pos1 < pos2 <= posb) 
	for i in posa:posb
		if i < pos1
			vd2 = id(leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		elseif i == pos1
			data2[i] = u_left
		elseif i < pos2
			vd2 = id(leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
		elseif i == pos2
			data2[i] = v_left
		else
			vd2 = id(leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
		end
		leftspace = space_r(data2[i])'
	end
	for j in lattice.k-1:-1:3
		posa, posb = band_boundary(lattice, j)
		pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
		@assert (posa <= pos1 < pos2 <= posb) 
		for i in posa:posb
			if i < pos1
				vd2 = id(leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
			elseif i == pos1
				data2[i] = u_middle
			elseif i < pos2
				vd2 = id(leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
			elseif i == pos2
				data2[i] = v_middle
			else
				vd2 = id(leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
			end
			leftspace = space_r(data2[i])'
		end
	end
	j = 2
	posa, posb = band_boundary(lattice, j)
	pos1, pos2 = index(lattice, j, conj=false, band=band), index(lattice, j, conj=true, band=band)
	@assert (posa <= pos1 < pos2 <= posb) 
	for i in posa:posb
		if i < pos1
			vd2 = id(leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
		elseif i == pos1
			data2[i] = u_right
		elseif i < pos2
			vd2 = id(leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
		elseif i == pos2
			data2[i] = v_right
		else
			vd2 = id(leftspace)
			data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]			
		end
		leftspace = space_r(data2[i])'
	end
	j = 1
	posa, posb = band_boundary(lattice, j)
	for pos in posa:posb
		vd2 = id(leftspace)
		data2[pos] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
	end
	return MPO(data2)
end

