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
function ti_mpotensor(corr::CorrelationMatrix, alg::ExponentialExpansionAlgorithm)
	ph = z2space()
	f = isomorphism(fuse(ph, ph), ph ⊗ ph)
	_JW = JW
	@tensor a[1,7;3,8] := σ₊[1,2,3,4] * _JW[5,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor abar[1,7;3,8] := _JW[2,4] * σ₋'[1,5,3,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor JW4[5;6] := _JW[1,2] * _JW[3,4] * f[5,1,3] * conj(f[6,2,4])
	m1 = GenericDecayTerm(a, abar, -corr.ηₖⱼ[2:end], middle = JW4)

	@tensor a[1,7;3,8] := σ₋'[1,2,3,4] * I2[5,6] * f[7,2,5] * conj(f[8,4,6])
	@tensor abar[1,7;3,8] := I2[2,4] * σ₊[1,5,3,6] * f[7,2,5] * conj(f[8,4,6])
	m2 = GenericDecayTerm(abar, a, corr.ηⱼₖ[2:end], middle = JW4)

	m1s = expand_decayterm(m1, alg=alg)
	m2s = expand_decayterm(m2, alg=alg)
	@tensor abar_a[7;8] := σ₊[1,2,3,4] * σ₋'[3,5,1,6] * f[7,2,5] * conj(f[8,4,6])

	coef = corr[1, 1]
	return SchurMPOTensor(-coef * abar_a, vcat(m1s, m2s))
end



