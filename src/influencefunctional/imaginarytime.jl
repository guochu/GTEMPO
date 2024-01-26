function partialinfluencefunctional(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector; band::Int=1)
	row = index(lattice, i, band=band, conj=true)
	col_pos = [index(lattice, j, band=band, conj=false) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end

"""
	influenceoperator(lattice::ImagGrassmannLattice1Order{<:A1A1B1B1}, corr2::ImagCorrelationFunction; band, algexpan)
"""
function influenceoperator(lattice::ImagGrassmannLattice1Order, corr2::ImagCorrelationFunction; band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	mpoj = ti_mpotensor(corr, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	mpo = MPO(h)
	_JW = JW
	return _fit_to_lattice(lattice, mpo, _JW, band) 
end

"""
	influenceoperatorexponential(lattice::ImagGrassmannLattice1Order{<:A1A1B1B1}, corr2::ImagCorrelationFunction, dt, alg; band, algexpan)
"""
function influenceoperatorexponential(lattice::ImagGrassmannLattice1Order, corr2::ImagCorrelationFunction, dt::Real, alg::FirstOrderStepper; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	mpoj = ti_mpotensor(corr, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	h2 = timeevompo(h, dt, alg)
	mpo = MPO(h2)
	_JW = I2
	return _fit_to_lattice(lattice, mpo, _JW, band) 
end
function influenceoperatorexponential(lattice::ImagGrassmannLattice1Order, corr2::ImagCorrelationFunction, dt::Real, alg::ComplexStepper; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	corr = corr2.data
	mpoj = ti_mpotensor(corr, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	h1, h2 = timeevompo(h, dt, alg)
	mpo1, mpo2 = MPO(h1), MPO(h2)
	_JW = I2
	return _fit_to_lattice(lattice, mpo1, _JW, band), _fit_to_lattice(lattice, mpo2, _JW, band) 
end

function differentialinfluencefunctional(lattice::ImagGrassmannLattice1Order{O}, corr::ImagCorrelationFunction, dt::Real, alg::TimeEvoMPOAlgorithm, algmult::DMRGAlgorithm; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion()) where O
	if !(LayoutStyle(lattice) isa TimeLocalLayout)
		lattice2 = similar(lattice, ordering = A1A1B1B1())
		mps = _differentialinfluencefunctional(lattice2, corr, dt, alg, algmult; band=band, algexpan=algexpan)
		_, mps2 = changeordering(O, lattice2, mps, trunc=algmult.trunc)
		return mps2
	else
		return _differentialinfluencefunctional(lattice, corr, dt, alg, algmult; band=band, algexpan=algexpan)
	end
end


function _differentialinfluencefunctional(lattice::ImagGrassmannLattice1Order, corr::ImagCorrelationFunction, dt::Real, alg::FirstOrderStepper, algmult::DMRGAlgorithm; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	mps = influenceoperatorexponential(lattice, corr, dt, alg, band=band, algexpan=algexpan) * vacuumstate(lattice)
	return mps
end
function _differentialinfluencefunctional(lattice::ImagGrassmannLattice1Order, corr::ImagCorrelationFunction, dt::Real, alg::ComplexStepper, algmult::DMRGAlgorithm; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
 	mpo1, mpo2 = influenceoperatorexponential(lattice, corr, dt, alg, band=band, algexpan=algexpan)
	mps1 = mpo1 * vacuumstate(lattice)
	mps2 = mpo2 * vacuumstate(lattice)
	return mult(mps1, mps2, algmult)
end

function _fit_to_lattice(lattice::ImagGrassmannLattice, mpo::MPO, _JW::MPSBondTensor, band::Int, trunc::TruncationScheme=DefaultMPOTruncation)
	@assert length(mpo) == 3
	(LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("only works lattice with TimeLocalLayout"))

	u_left, v_left = split_mpotensor(mpo[1], trunc)
	u_middle, v_middle = split_mpotensor(mpo[2], trunc)
	u_right, v_right = split_mpotensor(mpo[3], trunc)
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

