const _AllowedRealGrassmannOrdering = Union{A1Ā1a1ā1B1B̄1b1b̄1, A1Ā1B1B̄1a1ā1b1b̄1}

# full influence operator, only works for ordering A1Ā1a1ā1B1B̄1b1b̄1 and A1Ā1B1B̄1a1ā1b1b̄1
function influenceoperator(lattice::RealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction; 
							band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr, band)
	mpoj1, mpoj2, mpoj3, mpoj4 = ti_mpotensor(η⁺⁺, algexpan), ti_mpotensor(η⁺⁻, algexpan), ti_mpotensor(η⁻⁺, algexpan), ti_mpotensor(η⁻⁻, algexpan)
	# η⁺⁻ = CorrelationMatrix(η⁺⁻.ηₖⱼ, η⁺⁻.ηⱼₖ)
	# η⁻⁺ = CorrelationMatrix(η⁻⁺.ηₖⱼ, η⁻⁺.ηⱼₖ)
	h1, h2, h3, h4 = _get_mpo3(mpoj1), _get_mpo3(mpoj2), _get_mpo3(mpoj3), _get_mpo3(mpoj4)
	_JW = JW
	mpo1 = _fit_to_lattice(lattice, h1, _JW, band, :+, :+)
	# noticing the following two!!!
	mpo2 = _fit_to_lattice(lattice, h2, _JW, band, :-, :+) 
	mpo3 = _fit_to_lattice(lattice, h3, _JW, band, :+, :-) 
	mpo4 = _fit_to_lattice(lattice, h4, _JW, band, :-, :-) 
	return mpo1, mpo2, mpo3, mpo4
end

function influenceoperatorexponential(lattice::RealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::FirstOrderStepper; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr, band)
	mpoj1, mpoj2, mpoj3, mpoj4 = ti_mpotensor(η⁺⁺, algexpan), ti_mpotensor(η⁺⁻, algexpan), ti_mpotensor(η⁻⁺, algexpan), ti_mpotensor(η⁻⁻, algexpan)
	mpoj1, mpoj2, mpoj3, mpoj4 = timeevompo(mpoj1, dt, alg), timeevompo(mpoj2, dt, alg), timeevompo(mpoj3, dt, alg), timeevompo(mpoj4, dt, alg)
	h1, h2, h3, h4 = _get_mpo3(mpoj1), _get_mpo3(mpoj2), _get_mpo3(mpoj3), _get_mpo3(mpoj4)
	_JW = I2
	mpo1 = _fit_to_lattice(lattice, h1, _JW, band, :+, :+)
	# noticing the following two!!!
	mpo2 = _fit_to_lattice(lattice, h2, _JW, band, :-, :+) 
	mpo3 = _fit_to_lattice(lattice, h3, _JW, band, :+, :-) 
	mpo4 = _fit_to_lattice(lattice, h4, _JW, band, :-, :-) 
	return mpo1, mpo2, mpo3, mpo4
end
function influenceoperatorexponential(lattice::RealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::ComplexStepper; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr, band)
	mpoj1, mpoj2, mpoj3, mpoj4 = ti_mpotensor(η⁺⁺, algexpan), ti_mpotensor(η⁺⁻, algexpan), ti_mpotensor(η⁻⁺, algexpan), ti_mpotensor(η⁻⁻, algexpan)
	mpoj1a, mpoj1b = timeevompo(mpoj1, dt, alg)
	mpoj2a, mpoj2b = timeevompo(mpoj2, dt, alg)
	mpoj3a, mpoj3b = timeevompo(mpoj3, dt, alg)
	mpoj4a, mpoj4b = timeevompo(mpoj4, dt, alg)
	h1a, h1b = _get_mpo3(mpoj1a), _get_mpo3(mpoj1b)
	h2a, h2b = _get_mpo3(mpoj2a), _get_mpo3(mpoj2b)
	h3a, h3b = _get_mpo3(mpoj3a), _get_mpo3(mpoj3b)
	h4a, h4b = _get_mpo3(mpoj4a), _get_mpo3(mpoj4b)
	_JW = I2
	mpo1a, mpo1b = _fit_to_lattice(lattice, h1a, _JW, band, :+, :+), _fit_to_lattice(lattice, h1b, _JW, band, :+, :+)
	# noticing the following two!!!
	mpo2a, mpo2b = _fit_to_lattice(lattice, h2a, _JW, band, :-, :+), _fit_to_lattice(lattice, h2b, _JW, band, :-, :+) 
	mpo3a, mpo3b = _fit_to_lattice(lattice, h3a, _JW, band, :+, :-), _fit_to_lattice(lattice, h3a, _JW, band, :+, :-) 
	mpo4a, mpo4b = _fit_to_lattice(lattice, h4a, _JW, band, :-, :-), _fit_to_lattice(lattice, h4b, _JW, band, :-, :-) 
	return (mpo1a, mpo1b), (mpo2a, mpo2b), (mpo3a, mpo3b), (mpo4a, mpo4b)
end
function differentialinfluencefunctional(lattice::RealGrassmannLattice{O}, corr::RealCorrelationFunction, dt::Real, alg::TimeEvoMPOAlgorithm, algmult::DMRGAlgorithm;
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion()) where O
	if !(OrderingStyle(lattice) isa _AllowedRealGrassmannOrdering)
		lattice2 = similar(lattice, ordering = A1Ā1a1ā1B1B̄1b1b̄1())
		mps = _differentialinfluencefunctional(lattice2, corr, dt, alg, algmult; band=band, algexpan=algexpan)
		_, mps2 = changeordering(O, lattice2, mps, trunc=algmult.trunc)
		return mps2
	else
		return _differentialinfluencefunctional(lattice, corr, dt, alg, algmult; band=band, algexpan=algexpan)
	end
end

function _differentialinfluencefunctional(lattice::RealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::FirstOrderStepper, 
											algmult::DMRGAlgorithm; band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	h1, h2, h3, h4 = influenceoperatorexponential(lattice, corr, dt, alg, band=band, algexpan=algexpan)
	mps = h1 * vacuumstate(lattice)
	tmp = h2 * vacuumstate(lattice)
	mps = mult(mps, tmp, algmult)
	tmp = h3 * vacuumstate(lattice)
	mps = mult(mps, tmp, algmult)
	tmp = h4 * vacuumstate(lattice)
	mps = mult(mps, tmp, algmult)
	return mps
end
function _differentialinfluencefunctional(lattice::RealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::ComplexStepper, 
											algmult::DMRGAlgorithm; band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	(h1a, h1b), (h2a, h2b), (h3a, h3b), (h4a, h4b) = influenceoperatorexponential(lattice, corr, dt, alg, band=band, algexpan=algexpan)
	mps1 = h1a * vacuumstate(lattice)
	tmp = h1b * vacuumstate(lattice)
	mps = mult(mps1, tmp, algmult)

	mps1 = h2a * vacuumstate(lattice)
	tmp = h2b * vacuumstate(lattice)
	mps1 = mult(mps1, tmp, algmult)
	mps = mult(mps, mps1, algmult)

	mps1 = h3a * vacuumstate(lattice)
	tmp = h3b * vacuumstate(lattice)
	mps1 = mult(mps1, tmp, algmult)
	mps = mult(mps, mps1, algmult)

	mps1 = h4a * vacuumstate(lattice)
	tmp = h4b * vacuumstate(lattice)
	mps1 = mult(mps1, tmp, algmult)
	mps = mult(mps, mps1, algmult)

	return mps
end

function _get_mpo3(mpoj)
	# mpoj = ti_mpotensor(η, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	return MPO(h)
end
function _get_signed_corr(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, band::Int)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	if index(lattice, 1, conj=false, branch=:+, band=band) > index(lattice, 1, conj=true, branch=:+, band=band)
		η⁺⁺ = -transpose(η⁺⁺)
	end
	if index(lattice, 1, conj=false, branch=:-, band=band) > index(lattice, 1, conj=true, branch=:+, band=band)
		η⁺⁻ = -transpose(η⁺⁻)
	end
	if index(lattice, 1, conj=false, branch=:+, band=band) > index(lattice, 1, conj=true, branch=:-, band=band)
		η⁻⁺ = -transpose(η⁻⁺)
	end
	if index(lattice, 1, conj=false, branch=:-, band=band) > index(lattice, 1, conj=true, branch=:-, band=band)
		η⁻⁻ = -transpose(η⁻⁻)
	end
	return η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻
end

function _fit_to_lattice(lattice::RealGrassmannLattice, mpo::MPO, _JW::MPSBondTensor, band::Int, f1::Symbol, f2::Symbol, trunc::TruncationScheme=DefaultMPOTruncation)
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
	pos1, pos2 = index(lattice, j, conj=false, branch=f1, band=band), index(lattice, j, conj=true, branch=f2, band=band)
	if pos1 > pos2 # this sign has already been taken care of
		pos1, pos2 = pos2, pos1
	end
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
	for j in lattice.k-1:-1:2
		posa, posb = band_boundary(lattice, j)
		pos1, pos2 = index(lattice, j, conj=false, branch=f1, band=band), index(lattice, j, conj=true, branch=f2, band=band)
		if pos1 > pos2 # this sign has already been taken care of
			pos1, pos2 = pos2, pos1
		end
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
	j = 1
	posa, posb = band_boundary(lattice, j)
	pos1, pos2 = index(lattice, j, conj=false, branch=f1, band=band), index(lattice, j, conj=true, branch=f2, band=band)
	if pos1 > pos2 # this sign has already been taken care of
		pos1, pos2 = pos2, pos1
	end
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
	return MPO(data2)
end

