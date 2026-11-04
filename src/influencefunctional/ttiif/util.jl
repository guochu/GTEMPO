# build the translationally invariant MPOTensor of IF
"""
	ti_mpotensor(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm)

Return the 3-site unit cell of the translationally invariant MPO representing
`exp(Σᵢⱼ āᵢ Δᵢⱼ aⱼ)` on a chain whose site index increases into the past.
The block is read through `corr[i,j]`: `ηⱼₖ[i-j+1]` for i>j, `ηₖⱼ[j-i+1]`
for i<j and `ηⱼₖ[1]+ηₖⱼ[1]` at equal time.

The MPO lists the operators from left (future) to right (past), while the
influence-functional exponent is ordered `āᵢ Δᵢⱼ aⱼ`. A term that has `a` to
the left of `ā` therefore comes with a minus,
`āᵢ Δᵢⱼ aⱼ = -aᵢ Δⱼᵢ āⱼ`. This is why the i>j long-range term uses
`-ηₖⱼ[d]` (its `a` sits on the future leg) and the equal-time term uses
`-corr[1,1]`, while the i<j term uses `+ηⱼₖ[d]`.
"""
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# vertex-shifted bulk cell (2026-10, same convention as the QuAPI fix in
# models.jl / contourindex.jl): the strict node-vertex IF recombines the
# frozen cell-integral Toeplitz generators. Family/tail assignment (settled
# against the accessor kernel audit, GTEMPO/performance/firstorder/
# lag1_contact_audit.jl): the m1 family (a ⊗ JW-middles ⊗ ā) populates the
# upper lags and the m2 family (ā ⊗ middles ⊗ a) the lower lags of the
# same-branch blocks; the τ ring was reverted to the raw cell in 2026-10
# (the shift worsened its g(τ) discrete error constant ~30×) and the :τ
# branch is gone. Per block (kind = :+ / :-, the ηⱼₖ/ηₖⱼ roles are swapped
# between the two):
# - diag slot: ηₖⱼ[2] (the old upper family head);
# - upper lags: head-dropped [ηₖⱼ[3:M]] with NO far-continuation — the
#   continued last value is a spurious pair on a windowed turning-point leg
#   (removed by `_window_correction`, which only ever reads in-range slots)
#   and its repetition kink poisons the Prony fit of the whole family;
# - lower lags: the tail [ηⱼₖ[2:M-1]] delayed by one slot (fitted and
#   divided by the poles, so slot d ≥ 2 carries the tail value and slot 1
#   only the smooth backward extrapolation d1 = Σ cₖ of the raw fit), plus a
#   native `NearestNeighbourTerm` contact row with coeff E - d1: the delayed
#   channels and the contact row together equal the equal-time sum E exactly
#   at distance 1 and the tail unchanged beyond, with no dynamic-range
#   degeneracy (an earlier tiny-α channel carrying this correction
#   underflowed its coupling block and put ~1e14 on the diagonal).
# The cross blocks keep their values (their turning-point rows/cols are
# windowed to zero by the `_window_correction` factors).
# ---------------------------------------------------------------------------
function _vertex_shifted_parts(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm, kind::Symbol)
	E = corr[1, 1]
	M = length(corr.ηⱼₖ)
	M == 1 && return _ti_mpo_parts(corr, alg)
	ph = z2space()
	f = isomorphism(fuse(ph, ph), ph ⊗ ph)
	_JW = JW
	@tensor a1[1,7;3,8] := σ₊[1,2,3,4] * _JW[5,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor abar1[1,7;3,8] := _JW[2,4] * σ₋'[1,5,3,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor JW4[5;6] := _JW[1,2] * _JW[3,4] * f[5,1,3] * conj(f[6,2,4])
	@tensor a2[1,7;3,8] := σ₋'[1,2,3,4] * I2[5,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor abar2[1,7;3,8] := I2[2,4] * σ₊[1,5,3,6] * f[7,2,5] * conj(f[8,4,6])
	_delay(t) = ExponentialDecayTerm(t.a, t.b; middle=t.m, α=t.α,
									 coeff=iszero(t.coeff) ? t.coeff : t.coeff / t.α)
	if kind === :+
		# upper: head-dropped continuation, no far-continuation
		m1s = expand_decayterm(GenericDecayTerm(a1, abar1, -corr.ηₖⱼ[3:M], middle=JW4), alg=alg)
		# lower: tail delayed by one slot + native contact row
		Ltail = corr.ηⱼₖ[2:M-1]
		(length(Ltail) == 1) && (Ltail = vcat(Ltail, Ltail))   # Prony 下限
		m2s = expand_decayterm(GenericDecayTerm(abar2, a2, Ltail, middle=JW4), alg=alg)
		m2s = [_delay(t) for t in m2s]
		d1 = sum(t.coeff * t.α for t in m2s)   # 延迟族 slot-1 振幅 = Σ cₖ（除法前的拟合系数）
		m2s = vcat(m2s, [NearestNeighbourTerm(abar2, a2; coeff=E - d1)])
		return corr.ηₖⱼ[2], m1s, m2s
	else                                    # :- (mirror: the time-reversed --
		# block swaps the two decay families and the ηⱼₖ/ηₖⱼ roles)
		m2s = expand_decayterm(GenericDecayTerm(abar2, a2, corr.ηⱼₖ[3:M], middle=JW4), alg=alg)
		Utail = corr.ηₖⱼ[2:M-1]
		(length(Utail) == 1) && (Utail = vcat(Utail, Utail))
		m1s = expand_decayterm(GenericDecayTerm(a1, abar1, -Utail, middle=JW4), alg=alg)
		m1s = [_delay(t) for t in m1s]
		d1 = sum(t.coeff * t.α for t in m1s)   # 延迟族 slot-1 振幅 = Σ cₖ（除法前的拟合系数）
		# m1 家族的核贡献带算符序翻转（a 在未来腿：ā…a 核取 −），接触行同理
		m1s = vcat(m1s, [NearestNeighbourTerm(a1, abar1; coeff=-E - d1)])
		return corr.ηⱼₖ[2], m1s, m2s
	end
end

function _ti_mpotensor_shifted(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm, kind::Symbol)
	return _ti_mpotensor(_vertex_shifted_parts(corr, alg, kind)...)
end

# boundary variant of the shifted cell with individual slots switched off
# (zeroed coefficients keep the virtual-space structure, so the boundary cells
# glue into the same 3-cell `MPOHamiltonian` as the full middle cell)
function _ti_mpotensor_shifted_variant(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm, kind::Symbol;
									   onsite::Bool=true, m1::Bool=true, m2::Bool=true)
	coef, m1s, m2s = _vertex_shifted_parts(corr, alg, kind)
	return _ti_mpotensor(onsite ? coef : zero(coef), m1 ? m1s : _zerocoeff(m1s), m2 ? m2s : _zerocoeff(m2s))
end

# the three operator slots of the unit cell: the equal-node onsite block and
# the two decay-term families m1 = (a future, ā past) and m2 = (ā future, a
# past), expanded into exponential channels
function _ti_mpo_parts(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm)
	ph = z2space()
	f = isomorphism(fuse(ph, ph), ph ⊗ ph)
	_JW = JW
	@tensor a[1,7;3,8] := σ₊[1,2,3,4] * _JW[5,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor abar[1,7;3,8] := _JW[2,4] * σ₋'[1,5,3,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor JW4[5;6] := _JW[1,2] * _JW[3,4] * f[5,1,3] * conj(f[6,2,4])
	m1s = expand_decayterm(GenericDecayTerm(a, abar, -corr.ηₖⱼ[2:end], middle = JW4), alg=alg)

	@tensor a[1,7;3,8] := σ₋'[1,2,3,4] * I2[5,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor abar[1,7;3,8] := I2[2,4] * σ₊[1,5,3,6] * f[7,2,5] * conj(f[8,4,6])
	m2s = expand_decayterm(GenericDecayTerm(abar, a, corr.ηⱼₖ[2:end], middle = JW4), alg=alg)
	return corr[1, 1], m1s, m2s
end

function _ti_mpotensor(coef, m1s, m2s)
	ph = z2space()
	f = isomorphism(fuse(ph, ph), ph ⊗ ph)
	@tensor abar_a[7;8] := σ₊[1,2,3,4] * σ₋'[3,5,1,6] * f[7,2,5] * conj(f[8,4,6])
	return SchurMPOTensor(-coef * abar_a, vcat(m1s, m2s))
end

function ti_mpotensor(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm)
	coef, m1s, m2s = _ti_mpo_parts(corr, alg)
	return _ti_mpotensor(coef, m1s, m2s)
end

# boundary variant of the unit cell with individual slots switched off (the
# real-time vertex-shifted convention windows the turning-point legs: e.g. the
# ++ block drops its onsite and m1 slots at both contour ends). The dropped
# decay channels keep their operator structure with zero coefficients, so the
# variant has the same virtual spaces as the full cell and can be glued into
# the same 3-cell `MPOHamiltonian`; a cell with everything dropped is replaced
# by fillers in `_fit_to_lattice` instead.
_zerocoeff(ts) = [t isa ExponentialDecayTerm ?
                   ExponentialDecayTerm(t.a, t.b; middle=t.m, α=t.α, coeff=zero(t.coeff)) :
                   NearestNeighbourTerm(t.a, t.b; coeff=zero(t.coeff)) for t in ts]
function _ti_mpotensor_variant(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm;
							   onsite::Bool=true, m1::Bool=true, m2::Bool=true)
	coef, m1s, m2s = _ti_mpo_parts(corr, alg)
	return _ti_mpotensor(onsite ? coef : zero(coef), m1 ? m1s : _zerocoeff(m1s), m2 ? m2s : _zerocoeff(m2s))
end



