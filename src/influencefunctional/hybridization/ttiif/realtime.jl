const _AllowedRealGrassmannOrdering = Union{A1Ā1a1ā1B1B̄1b1b̄1, A1Ā1B1B̄1a1ā1b1b̄1, a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2}

# full influence operator, works for the time-local adjacent-conjugation orderings
# A1Ā1a1ā1B1B̄1b1b̄1 and A1Ā1B1B̄1a1ā1b1b̄1 (time descending),
# a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2 (time ascending)
function influenceoperators(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction; 
							band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr, band)
	mpoj1, mpoj2, mpoj3, mpoj4 = ti_mpotensor(η⁺⁺, alg_expan), ti_mpotensor(η⁺⁻, alg_expan), ti_mpotensor(η⁻⁺, alg_expan), ti_mpotensor(η⁻⁻, alg_expan)
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

function influenceoperatorsteppers(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::FirstOrderStepper; 
										band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr, band)
	mpoj1, mpoj2, mpoj3, mpoj4 = ti_mpotensor(η⁺⁺, alg_expan), ti_mpotensor(η⁺⁻, alg_expan), ti_mpotensor(η⁻⁺, alg_expan), ti_mpotensor(η⁻⁻, alg_expan)
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
function influenceoperatorsteppers(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::ComplexStepper; 
										band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr, band)
	mpoj1, mpoj2, mpoj3, mpoj4 = ti_mpotensor(η⁺⁺, alg_expan), ti_mpotensor(η⁺⁻, alg_expan), ti_mpotensor(η⁻⁺, alg_expan), ti_mpotensor(η⁻⁻, alg_expan)
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
function influenceoperatorstepper(lattice::AbstractRealGrassmannLattice{O}, corr::RealCorrelationFunction, dt::Real, alg::TimeEvoMPOAlgorithm, alg_mult::DMRGAlgorithm;
										band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan()) where O
	if !(OrderingStyle(lattice) isa _AllowedRealGrassmannOrdering)
		lattice2 = similar(lattice, ordering = A1Ā1a1ā1B1B̄1b1b̄1())
		mps = _influenceoperatorstepper(lattice2, corr, dt, alg, alg_mult; band=band, alg_expan=alg_expan)
		_, mps2 = changeordering(O, lattice2, mps, trunc=alg_mult.trunc)
		return mps2
	else
		return _influenceoperatorstepper(lattice, corr, dt, alg, alg_mult; band=band, alg_expan=alg_expan)
	end
end

function _influenceoperatorstepper(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::FirstOrderStepper, 
											alg_mult::DMRGAlgorithm; band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan())
	h1, h2, h3, h4 = influenceoperatorsteppers(lattice, corr, dt, alg, band=band, alg_expan=alg_expan)
	mps = h1 * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	tmp = h2 * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	mps = mult(mps, tmp, alg_mult)[1]
	tmp = h3 * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	mps = mult(mps, tmp, alg_mult)[1]
	tmp = h4 * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	mps = mult(mps, tmp, alg_mult)[1]
	return mps
end
function _influenceoperatorstepper(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::ComplexStepper, 
											alg_mult::DMRGAlgorithm; band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan())
	(h1a, h1b), (h2a, h2b), (h3a, h3b), (h4a, h4b) = influenceoperatorsteppers(lattice, corr, dt, alg, band=band, alg_expan=alg_expan)
	mps1 = h1a * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	tmp = h1b * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	mps = mult(mps1, tmp, alg_mult)[1]

	mps1 = h2a * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	tmp = h2b * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	mps1 = mult(mps1, tmp, alg_mult)[1]
	mps = mult(mps, mps1, alg_mult)[1]

	mps1 = h3a * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	tmp = h3b * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	mps1 = mult(mps1, tmp, alg_mult)[1]
	mps = mult(mps, mps1, alg_mult)[1]

	mps1 = h4a * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	tmp = h4b * vacuumstate(promote_type(scalartype(lattice), scalartype(corr)), lattice)
	mps1 = mult(mps1, tmp, alg_mult)[1]
	mps = mult(mps, mps1, alg_mult)[1]

	return mps
end

function _get_mpo3(mpoj)
	# mpoj = ti_mpotensor(η, alg_expan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	return MPO(h)
end
function _get_signed_corr(lattice::AbstractRealGrassmannLattice, corr::RealCorrelationFunction, band::Int)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	if lattice[GrassmannContourIndex(1, :+, band, false)] > lattice[GrassmannContourIndex(1, :+, band, true)]
		η⁺⁺ = -transpose(η⁺⁺)
	end
	if lattice[GrassmannContourIndex(1, :-, band, false)] > lattice[GrassmannContourIndex(1, :+, band, true)]
		η⁺⁻ = -transpose(η⁺⁻)
	end
	if lattice[GrassmannContourIndex(1, :+, band, false)] > lattice[GrassmannContourIndex(1, :-, band, true)]
		η⁻⁺ = -transpose(η⁻⁺)
	end
	if lattice[GrassmannContourIndex(1, :-, band, false)] > lattice[GrassmannContourIndex(1, :-, band, true)]
		η⁻⁻ = -transpose(η⁻⁻)
	end
	return η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻
end

function _fit_to_lattice(lattice::AbstractRealGrassmannLattice, mpo::MPO, _JW::MPSBondTensor, band::Int, f1::Symbol, f2::Symbol, trunc::TruncationScheme=Defaults.alg_trunc_if())
	@assert length(mpo) == 3
	(LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("only works lattice with TimeLocalLayout"))

	u_left, v_left = split_mpotensor(mpo[1], trunc)
	u_middle, v_middle = split_mpotensor(mpo[2], trunc)
	u_right, v_right = split_mpotensor(mpo[3], trunc)
	L = length(lattice)
	data2 = Vector{typeof(u_middle)}(undef, L)
	leftspace = oneunit(z2space())

	# fill the time blocks in physical site order: on a time-descending
	# lattice the first block is j=k, while on a time-ascending lattice it
	# is j=1. The 3-site bath MPO starts with (u_left, v_left) at the left
	# boundary and ends with (u_right, v_right) at the far end.
	ascending = RealTimeOrderingStyle(lattice) isa TimeAscending
	jfirst = ascending ? 1 : lattice.k
	jlast  = ascending ? lattice.k : 1
	jmids  = ascending ? (2:lattice.k-1) : (lattice.k-1:-1:2)

	function placeblock!(j, u, v, leftspace)
		posa, posb = band_boundary(lattice, j)
		pos1, pos2 = lattice[GrassmannContourIndex(j, f1, band, false)], lattice[GrassmannContourIndex(j, f2, band, true)]
		if pos1 > pos2 # this sign has already been taken care of
			pos1, pos2 = pos2, pos1
		end
		for i in posa:posb
			if i < pos1
				vd2 = id(leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
			elseif i == pos1
				data2[i] = u
			elseif i < pos2
				vd2 = id(leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * _JW[3,4]
			elseif i == pos2
				data2[i] = v
			else
				vd2 = id(leftspace)
				data2[i] = @tensor tmp[1,3;2,4] := vd2[1,2] * I2[3,4]
			end
			leftspace = space_r(data2[i])'
		end
		return leftspace
	end

	leftspace = placeblock!(jfirst, u_left, v_left, leftspace)
	for j in jmids
		leftspace = placeblock!(j, u_middle, v_middle, leftspace)
	end
	leftspace = placeblock!(jlast, u_right, v_right, leftspace)
	return MPO(data2)
end

