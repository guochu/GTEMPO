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
# frozen cell-integral Toeplitz generators. Slot mapping of `ti_mpotensor`
# (established numerically): the m1 family (a ⊗ JW-middles ⊗ ā) populates the
# node-lower lags (ā later) and the m2 family (ā ⊗ middles ⊗ a) the
# node-upper lags (a later). The shifted bulk matrices are
# - `:τ` ring (n nodes): diag ηₖⱼ[2], lower lags [E, ηⱼₖ[2:n-1]] with the
#   equal-time sum E = corr[1,1], upper lag d = -lower[n-d] (antiperiodic
#   wrap);
# - `:+` real-time same-branch (column shift +1): diag ηₖⱼ[2], lower
#   [E, ηⱼₖ[2:M-1]], upper [ηₖⱼ[3:M], cont]; the `:-` block is the mirror
#   (diag ηⱼₖ[2], lower [ηⱼₖ[3:M], cont], upper [E, ηₖⱼ[2:M-1]]). The
#   cross blocks keep their values (their turning-point rows/cols are
#   windowed to zero by `_window_drops` instead).
# The nearest-neighbour E pair cannot ride the smooth Prony tails (it would
# break the exponential fit), so it is emitted as an explicit α = 0 channel;
# the farthest lag has no stored value and is continued with the last one
# (the spurious farthest pair sits on a windowed turning-point leg and is
# dropped by `_window_drops`).
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
	if kind === :τ || kind === :+
		# lower lags: the tail ηⱼₖ[2:M-1] is Prony-expanded as it is, and the
		# equal-time sum E at distance 1 is delivered by an extra tiny-α
		# channel that corrects the smooth backward extrapolation
		# Σcₖαₖ of the tail to E (value c·α^d: the end operator of the Schur
		# assembly carries the factor α, so α = 0 would kill the channel)
		# NOTE(2026-10-30): the tiny-α equal-time channel below is DISABLED —
		# it degenerates in the Schur/MPO assembly (the end operator carries α,
		# so at α = 1e-3 the coupling block underflows while the coeff
		# (E+d1)/αε ~ 1e14 lands on the diagonal), which made XTRGIF return
		# NaN and TDVPIF hit Infs. Without the channel the d=1 value is the
		# smooth Prony extrapolation d1 instead of the equal-time sum E —
		# the known ~1-5% short-chain inaccuracy, to be fixed with the
		# short-time Prony redesign.
		Ltail = corr.ηⱼₖ[2:M-1]
		(length(Ltail) == 1) && (Ltail = vcat(Ltail, Ltail))   # Prony 下限
		m1s = expand_decayterm(GenericDecayTerm(a1, abar1, -Ltail, middle = JW4), alg=alg)
		Uvec = vcat(corr.ηₖⱼ[3:M], corr.ηₖⱼ[M])
		(length(Uvec) == 1) && (Uvec = vcat(Uvec, Uvec))
		m2s = expand_decayterm(GenericDecayTerm(abar2, a2, Uvec, middle = JW4), alg=alg)
		return corr.ηₖⱼ[2], m1s, m2s
	else                                    # :- (mirror: the time-reversed --
		# block swaps the two decay families — m1 populates the upper lags and
		# m2 the lower lags — established numerically)
		# (tiny-α equal-time channel disabled — see the note in the :τ/:+ branch)
		Utail = corr.ηₖⱼ[2:M-1]
		(length(Utail) == 1) && (Utail = vcat(Utail, Utail))
		m1s = expand_decayterm(GenericDecayTerm(a1, abar1, -Utail, middle = JW4), alg=alg)
		Ltail = vcat(corr.ηⱼₖ[3:M], corr.ηⱼₖ[M])
		(length(Ltail) == 1) && (Ltail = vcat(Ltail, Ltail))
		m2s = expand_decayterm(GenericDecayTerm(abar2, a2, Ltail, middle = JW4), alg=alg)
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
_zerocoeff(ts) = [ExponentialDecayTerm(t.a, t.b; middle=t.m, α=t.α, coeff=zero(t.coeff)) for t in ts]
function _ti_mpotensor_variant(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm;
							   onsite::Bool=true, m1::Bool=true, m2::Bool=true)
	coef, m1s, m2s = _ti_mpo_parts(corr, alg)
	return _ti_mpotensor(onsite ? coef : zero(coef), m1 ? m1s : _zerocoeff(m1s), m2 ? m2s : _zerocoeff(m2s))
end



