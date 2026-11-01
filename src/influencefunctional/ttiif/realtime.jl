const _AllowedRealGrassmannOrdering = Union{A1Ā1a1ā1B1B̄1b1b̄1, A1Ā1B1B̄1a1ā1b1b̄1, a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2}

# full influence operator, works for the time-local adjacent-conjugation orderings
# A1Ā1a1ā1B1B̄1b1b̄1 and A1Ā1B1B̄1a1ā1b1b̄1 (time descending),
# a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2 (time ascending)
function influenceoperators(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction;
							band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan())
	ηs = _get_signed_corr(lattice, corr, band)
	# block identities (row ā branch, col a branch); the fit labels below are
	# the swapped (col, row) pairs of the transposed cross blocks
	blocks = ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
	_JW = JW
	return map(enumerate(blocks)) do (n, (b1, b2))
		mid = (b1 === b2) ? _ti_mpotensor_shifted(ηs[n], alg_expan, b1) : ti_mpotensor(ηs[n], alg_expan)
		h = MPO(MPOHamiltonian([mid, mid, mid]))
		_fit_to_lattice(lattice, h, _JW, band, b2, b1)
	end
end

function influenceoperatorsteppers(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::FirstOrderStepper;
										band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan())
	ηs = _get_signed_corr(lattice, corr, band)
	blocks = ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
	_JW = I2
	return map(enumerate(blocks)) do (n, (b1, b2))
		mid = (b1 === b2) ? _ti_mpotensor_shifted(ηs[n], alg_expan, b1) : ti_mpotensor(ηs[n], alg_expan)
		h = MPO(MPOHamiltonian([timeevompo(mid, dt, alg) for _ in 1:3]))
		_fit_to_lattice(lattice, h, _JW, band, b2, b1)
	end
end
function influenceoperatorsteppers(lattice::AbstractRealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, dt::Real, alg::ComplexStepper;
										band::Int=1, alg_expan::ExponentialExpansionAlgorithm=Defaults.alg_expan())
	ηs = _get_signed_corr(lattice, corr, band)
	blocks = ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
	_JW = I2
	return map(enumerate(blocks)) do (n, (b1, b2))
		mid = (b1 === b2) ? _ti_mpotensor_shifted(ηs[n], alg_expan, b1) : ti_mpotensor(ηs[n], alg_expan)
		c1, c2 = timeevompo(mid, dt, alg)
		h1, h2 = MPO(MPOHamiltonian([c1, c1, c1])), MPO(MPOHamiltonian([c2, c2, c2]))
		_fit_to_lattice(lattice, h1, _JW, band, b2, b1),
		_fit_to_lattice(lattice, h2, _JW, band, b2, b1)
	end
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

function _fit_to_lattice(lattice::AbstractRealGrassmannLattice, mpo::MPO, _JW::MPSBondTensor, band::Int, f1::Symbol, f2::Symbol,
						 trunc::TruncationScheme=Defaults.alg_trunc_if(); dropleft::Bool=false, dropright::Bool=false)
	@assert length(mpo) == 3
	(LayoutStyle(lattice) isa TimeLocalLayout) || throw(ArgumentError("only works lattice with TimeLocalLayout"))

	u_left, v_left = split_mpotensor(mpo[1], notrunc())
	u_middle, v_middle = split_mpotensor(mpo[2], notrunc())
	u_right, v_right = split_mpotensor(mpo[3], notrunc())
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

function _get_signed_corr(lattice::AbstractRealGrassmannLattice, corr::RealCorrelationFunction, band::Int)
	# The MPO built by `ti_mpotensor` writes the operators from left (future)
	# to right (past). If the column leg `a` of a block lies to the right of
	# its row leg `ā` at equal time, the block is stored in the opposite
	# (i,j) orientation: transpose it, and pick up a minus from
	# `āᵢ aⱼ = -āⱼ aᵢ`. The transposed block is then fitted under the swapped
	# branch labels in `influenceoperators`.
	# Vertex-shifted convention (2026-10): the same-branch blocks are built
	# from the shifted bulk generators at the parts level (see
	# `_ti_mpotensor_shifted`); the cross blocks keep their values.
	function signed(b1::Symbol, b2::Symbol)
		η = _physical_branch(corr, b1, b2)
		swap = lattice[GrassmannContourIndex(1, b2, band, false)] >
			   lattice[GrassmannContourIndex(1, b1, band, true)]
		return swap ? -transpose(η) : η
	end
	return signed(:+, :+), signed(:+, :-), signed(:-, :+), signed(:-, :-)
end

# ---------------------------------------------------------------------------
# boundary window correction of the vertex-shifted real-time convention: the
# strict node-vertex IF windows the turning-point legs (row ā⁺ / col a⁻ at
# t = 0 and col a⁺ / row ā⁻ at t = T) to zero, while the TTI bulk MPO built
# from the shifted generators carries the Toeplitz continuation there. Since
# all factors are exponentials of (ξ, η)-bilinears, they commute and the exact
# windowed IF = (TTI bulk MPO) × exp(-R), with R the bulk entries on the
# windowed legs. The correction factors are rank-1 partial IFs built with
# `partialmpo` and multiplied in by the consumer entry points (XTRGIF /
# ExactTTIIF / TDVPIF).
# ---------------------------------------------------------------------------
# correction factors exp(-R) of one branch block: R carries the bulk (Toeplitz
# continuation) entries on the windowed legs, so that (TTI bulk MPO) × exp(-R)
# = the strict windowed IF exactly (all factors are exponentials of commuting
# (ξ, η)-bilinears). The row / column factors are built exactly like the
# partial IF of `partialif_hybrid` (analytic bond-2 MPOs).
function _rowfactor(lattice, rowidx::GrassmannContourIndex, colidxs, coefs)
	row = lattice[rowidx]
	col_pos = [lattice[c] for c in colidxs]
	return partialmpo(row, col_pos, coefs) *
		   vacuumstate(promote_type(scalartype(lattice), eltype(coefs)), lattice)
end
function _colfactor(lattice, colidx::GrassmannContourIndex, rowidxs, coefs)
	col = lattice[colidx]
	row_pos = [lattice[r] for r in rowidxs]
	return partialmpo(col, row_pos, -coefs) *
		   vacuumstate(promote_type(scalartype(lattice), eltype(coefs)), lattice)
end

function _window_correction(lattice::AbstractRealGrassmannLattice, corr::RealCorrelationFunction,
							band::Int, b1::Symbol, b2::Symbol)
	k = lattice.k
	η = _physical_branch(corr, b1, b2)
	M = length(η.ηⱼₖ)
	E = η[1, 1]
	ridxs(b) = [GrassmannContourIndex(i, b, band, true) for i in 1:k]
	cidxs(b) = [GrassmannContourIndex(j, b, band, false) for j in 1:k]
	if b1 === :+ && b2 === :+
		U(d) = d + 2 <= M ? η.ηₖⱼ[d+2] : η.ηₖⱼ[M]
		dg = η.ηₖⱼ[2]
		# row ā⁺₁ / column a⁺_k removal. `partialmpo` attaches the
		# coefficients to the [listed legs × single leg] pairs (established
		# numerically), so the row sweep is realized as a column-style factor
		# on the a⁺_k leg and vice versa
		Rrow = [dg; [U(j-1) for j in 2:k]]                 # B[1, j], ā⁺₁ row
		Rcol = [i == k ? dg : U(k-i) for i in 1:k]         # B[i, k], a⁺_k col
		Rcor = U(k-1)
		# 行扫（ā⁺₁ 行，右端腿在序外）：逐对因子；列扫（a⁺_k 列，左端在序内）：
		# partialmpo 列式因子
		fs = GrassmannMPS[]
		for j in 1:k
			push!(fs, partialmpo(lattice[GrassmannContourIndex(j, :+, band, false)],
								 [lattice[GrassmannContourIndex(1, :+, band, true)]], [-Rrow[j]]) *
					   vacuumstate(scalartype(corr), lattice))
		end
		push!(fs, _colfactor(lattice, GrassmannContourIndex(k, :+, band, false), ridxs(:+), -Rcol))
		push!(fs, partialmpo(lattice[GrassmannContourIndex(k, :+, band, false)],
							 [lattice[GrassmannContourIndex(1, :+, band, true)]], [Rcor]) *
				   vacuumstate(scalartype(corr), lattice))
		return fs
	elseif b1 === :- && b2 === :-
		L(d) = d + 2 <= M ? η.ηⱼₖ[d+2] : η.ηⱼₖ[M]
		dg = η.ηⱼₖ[2]
		Rrow = [j == k ? dg : L(k-j) for j in 1:k]         # B[k, j], ā⁻_k row
		Rcol = [dg; [L(i-1) for i in 2:k]]                 # B[i, 1], a⁻_1 col
		Rcor = L(k-1)
		frow = _rowfactor(lattice, GrassmannContourIndex(k, :-, band, true), cidxs(:-), -Rrow)
		fcol = _colfactor(lattice, GrassmannContourIndex(1, :-, band, false), ridxs(:-), -Rcol)
		fcor = _rowfactor(lattice, GrassmannContourIndex(k, :-, band, true),
						  [GrassmannContourIndex(1, :-, band, false)], [Rcor])
		return [frow, fcol, fcor]
	elseif b1 === :+ && b2 === :-
		Rrow = [E; [η.ηₖⱼ[j] for j in 2:k]]                # B[1, j], ā⁺₁ row
		Rcol = [E; [η.ηⱼₖ[i] for i in 2:k]]                # B[i, 1], a⁻_1 col
		frow = _rowfactor(lattice, GrassmannContourIndex(1, :+, band, true), cidxs(:-), -Rrow)
		fcol = _colfactor(lattice, GrassmannContourIndex(1, :-, band, false), ridxs(:+), -Rcol)
		fcor = _rowfactor(lattice, GrassmannContourIndex(1, :+, band, true),
						  [GrassmannContourIndex(1, :-, band, false)], [E])
		return [frow, fcol, fcor]
	else                                                    # (-, +): window at node k
		Rrow = [j == k ? E : η.ηⱼₖ[k-j+1] for j in 1:k]   # B[k, j], ā⁻_k row
		Rcol = [i == k ? E : η.ηₖⱼ[k-i+1] for i in 1:k]   # B[i, k], a⁺_k col
		frow = _rowfactor(lattice, GrassmannContourIndex(k, :-, band, true), cidxs(:+), -Rrow)
		fcol = _colfactor(lattice, GrassmannContourIndex(k, :+, band, false), ridxs(:-), -Rcol)
		fcor = _rowfactor(lattice, GrassmannContourIndex(k, :-, band, true),
						  [GrassmannContourIndex(k, :+, band, false)], [E])
		return [frow, fcol, fcor]
	end
end

# all four blocks' correction factors (multiplied in by the consumer entry
# points after the bulk TTI construction; the factors mutually commute)
function _window_corrections(lattice::AbstractRealGrassmannLattice, corr::RealCorrelationFunction; band::Int=1)
	blocks = ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
	return reduce(vcat, [_window_correction(lattice, corr, band, b1, b2) for (b1, b2) in blocks])
end
