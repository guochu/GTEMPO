function partialinfluencefunctional(lattice::RealGrassmannLattice, i::Int, cols::AbstractVector; fi::Bool, fj::Bool, band::Int=1)
	row = index(lattice, i, band=band, conj=true, forward=fi)
	col_pos = [index(lattice, j, band=band, conj=false, forward=fj) for j in length(cols):-1:1]
	mpo = partialmpo(row, col_pos, reverse(cols))
	return mpo * vacuumstate(lattice)
end
function partialinfluencefunctional(lattice::RealGrassmannLattice, i::Int, cols_f::AbstractVector, cols_b::AbstractVector; fi::Bool, band::Int=1)
	row = index(lattice, i, band=band, conj=true, forward=fi)
	cols = eltype(cols_f)[]
	col_pos = Int[]
	for j in length(cols_f):-1:1
		pos1, pos2 = index(lattice, j, band=band, conj=false, forward=true), index(lattice, j, band=band, conj=false, forward=false)
		push!(col_pos, pos1)
		push!(col_pos, pos2)
		push!(cols, cols_f[j])
		push!(cols, cols_b[j])
	end
	mpo = partialmpo(row, col_pos, cols)
	return mpo * vacuumstate(lattice)
end
function partialinfluencefunctional(lattice::RealGrassmannLattice, rows::AbstractVector, j::Int; fi::Bool, fj::Bool, band::Int=1)
	col = index(lattice, j, band=band, conj=false, forward=fj)
	row_pos = [index(lattice, i, band=band, conj=true, forward=fi) for i in length(rows):-1:1]
	mpo = partialmpo(col, row_pos, -reverse(rows))
	return mpo * vacuumstate(lattice)
end
function partialinfluencefunctional(lattice::RealGrassmannLattice, rows_f::AbstractVector, rows_b::AbstractVector, j::Int; fj::Bool, band::Int=1)
	col = index(lattice, j, band=band, conj=false, forward=fj)
	rows = eltype(rows_f)[]
	row_pos = Int[]
	for i in length(rows_f):-1:1
		pos1, pos2 = index(lattice, i, band=band, conj=true, forward=true), index(lattice, i, band=band, conj=true, forward=false)
		push!(row_pos, pos1)
		push!(row_pos, pos2)
		push!(rows, -rows_f[i])
		push!(rows, -rows_b[i])
	end
	mpo = partialmpo(col, row_pos, rows)
	return mpo * vacuumstate(lattice)
end

# full influence operator, only works for ordering A1A1a1a1B1B1b1b1
function influenceoperator(lattice::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, corr::RealCorrelationFunction; band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr, band)
	mpoj1, mpoj2, mpoj3, mpoj4 = ti_mpotensor(η⁺⁺, algexpan), ti_mpotensor(η⁺⁻, algexpan), ti_mpotensor(η⁻⁺, algexpan), ti_mpotensor(η⁻⁻, algexpan)
	# η⁺⁻ = CorrelationMatrix(η⁺⁻.ηₖⱼ, η⁺⁻.ηⱼₖ)
	# η⁻⁺ = CorrelationMatrix(η⁻⁺.ηₖⱼ, η⁻⁺.ηⱼₖ)
	h1, h2, h3, h4 = _get_mpo3(mpoj1), _get_mpo3(mpoj2), _get_mpo3(mpoj3), _get_mpo3(mpoj4)
	_JW = JW
	mpo1 = _fit_to_lattice(lattice, h1, _JW, band, true, true)
	# noticing the following two!!!
	mpo2 = _fit_to_lattice(lattice, h2, _JW, band, false, true) 
	mpo3 = _fit_to_lattice(lattice, h3, _JW, band, true, false) 
	mpo4 = _fit_to_lattice(lattice, h4, _JW, band, false, false) 
	return mpo1, mpo2, mpo3, mpo4
end

function influenceoperatorexponential(lattice::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, corr::RealCorrelationFunction, dt::Real, alg::FirstOrderStepper; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = _get_signed_corr(lattice, corr, band)
	mpoj1, mpoj2, mpoj3, mpoj4 = ti_mpotensor(η⁺⁺, algexpan), ti_mpotensor(η⁺⁻, algexpan), ti_mpotensor(η⁻⁺, algexpan), ti_mpotensor(η⁻⁻, algexpan)
	mpoj1, mpoj2, mpoj3, mpoj4 = timeevompo(mpoj1, dt, alg), timeevompo(mpoj2, dt, alg), timeevompo(mpoj3, dt, alg), timeevompo(mpoj4, dt, alg)
	h1, h2, h3, h4 = _get_mpo3(mpoj1), _get_mpo3(mpoj2), _get_mpo3(mpoj3), _get_mpo3(mpoj4)
	_JW = I2
	mpo1 = _fit_to_lattice(lattice, h1, _JW, band, true, true)
	# noticing the following two!!!
	mpo2 = _fit_to_lattice(lattice, h2, _JW, band, false, true) 
	mpo3 = _fit_to_lattice(lattice, h3, _JW, band, true, false) 
	mpo4 = _fit_to_lattice(lattice, h4, _JW, band, false, false) 
	return mpo1, mpo2, mpo3, mpo4
end
function influenceoperatorexponential(lattice::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, corr::RealCorrelationFunction, dt::Real, alg::ComplexStepper; 
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
	mpo1a, mpo1b = _fit_to_lattice(lattice, h1a, _JW, band, true, true), _fit_to_lattice(lattice, h1b, _JW, band, true, true)
	# noticing the following two!!!
	mpo2a, mpo2b = _fit_to_lattice(lattice, h2a, _JW, band, false, true), _fit_to_lattice(lattice, h2b, _JW, band, false, true) 
	mpo3a, mpo3b = _fit_to_lattice(lattice, h3a, _JW, band, true, false), _fit_to_lattice(lattice, h3a, _JW, band, true, false) 
	mpo4a, mpo4b = _fit_to_lattice(lattice, h4a, _JW, band, false, false), _fit_to_lattice(lattice, h4b, _JW, band, false, false) 
	return (mpo1a, mpo1b), (mpo2a, mpo2b), (mpo3a, mpo3b), (mpo4a, mpo4b)
end
function differentialinfluencefunctional(lattice::RealGrassmannLattice{O}, corr::RealCorrelationFunction, dt::Real, alg::TimeEvoMPOAlgorithm;
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion(), trunc::TruncationScheme=DefaultMPOTruncation) where O
	if !(OrderingStyle(lattice) isa A1A1a1a1B1B1b1b1)
		lattice2 = similar(lattice, ordering = A1A1a1a1B1B1b1b1())
		mps = _differentialinfluencefunctional(lattice2, corr, dt, alg; band=band, algexpan=algexpan, trunc=trunc)
		_, mps2 = changeordering(O, lattice2, mps, trunc=trunc)
		return mps2
	else
		return _differentialinfluencefunctional(lattice, corr, dt, alg; band=band, algexpan=algexpan, trunc=trunc)
	end
end

function _differentialinfluencefunctional(lattice::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, corr::RealCorrelationFunction, dt::Real, alg::FirstOrderStepper; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion(), trunc::TruncationScheme=DefaultMPOTruncation)
	h1, h2, h3, h4 = influenceoperatorexponential(lattice, corr, dt, alg, band=band, algexpan=algexpan)
	mps = h1 * vacuumstate(lattice)
	tmp = h2 * vacuumstate(lattice)
	orth = Orthogonalize(trunc=trunc, normalize=false)
	mps = mult!(mps, tmp, trunc=trunc)
	tmp = h3 * vacuumstate(lattice)
	mps = mult!(mps, tmp, trunc=trunc)
	tmp = h4 * vacuumstate(lattice)
	mps = mult!(mps, tmp, trunc=trunc)
	return mps
end
function _differentialinfluencefunctional(lattice::RealGrassmannLattice{<:A1A1a1a1B1B1b1b1}, corr::RealCorrelationFunction, dt::Real, alg::ComplexStepper; 
										band::Int=1, algexpan::ExponentialExpansionAlgorithm=PronyExpansion(), trunc::TruncationScheme=DefaultMPOTruncation)
	(h1a, h1b), (h2a, h2b), (h3a, h3b), (h4a, h4b) = influenceoperatorexponential(lattice, corr, dt, alg, band=band, algexpan=algexpan)
	orth = Orthogonalize(trunc=trunc, normalize=false)
	mps = h1a * vacuumstate(lattice)
	canonicalize!(mps, alg=orth)
	mps = h1b * mps
	canonicalize!(mps, alg=orth)

	tmp = h2a * vacuumstate(lattice)
	canonicalize!(tmp, alg=orth)
	tmp = h2b * tmp
	canonicalize!(tmp, alg=orth)
	mps = mult!(mps, tmp, trunc=trunc)

	tmp = h3a * vacuumstate(lattice)
	canonicalize!(tmp, alg=orth)
	tmp = h3b * tmp
	canonicalize!(tmp, alg=orth)
	mps = mult!(mps, tmp, trunc=trunc)

	tmp = h4a * vacuumstate(lattice)
	canonicalize!(tmp, alg=orth)
	tmp = h4b * tmp
	canonicalize!(tmp, alg=orth)
	mps = mult!(mps, tmp, trunc=trunc)

	return mps
end

function _get_mpo3(mpoj)
	# mpoj = ti_mpotensor(η, algexpan)
	h = MPOHamiltonian([mpoj, mpoj, mpoj])
	return MPO(h)
end
function _get_signed_corr(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, band::Int)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	if index(lattice, 1, conj=false, forward=true, band=band) > index(lattice, 1, conj=true, forward=true, band=band)
		η⁺⁺ = -transpose(η⁺⁺)
	end
	if index(lattice, 1, conj=false, forward=false, band=band) > index(lattice, 1, conj=true, forward=true, band=band)
		η⁺⁻ = -transpose(η⁺⁻)
	end
	if index(lattice, 1, conj=false, forward=true, band=band) > index(lattice, 1, conj=true, forward=false, band=band)
		η⁻⁺ = -transpose(η⁻⁺)
	end
	if index(lattice, 1, conj=false, forward=false, band=band) > index(lattice, 1, conj=true, forward=false, band=band)
		η⁻⁻ = -transpose(η⁻⁻)
	end
	return η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻
end

function _fit_to_lattice(lattice::RealGrassmannLattice, mpo::MPO, _JW::MPSBondTensor, band::Int, f1::Bool, f2::Bool, trunc::TruncationScheme=DefaultMPOTruncation)
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
	pos1, pos2 = index(lattice, j, conj=false, forward=f1, band=band), index(lattice, j, conj=true, forward=f2, band=band)
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
		pos1, pos2 = index(lattice, j, conj=false, forward=f1, band=band), index(lattice, j, conj=true, forward=f2, band=band)
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
	pos1, pos2 = index(lattice, j, conj=false, forward=f1, band=band), index(lattice, j, conj=true, forward=f2, band=band)
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

