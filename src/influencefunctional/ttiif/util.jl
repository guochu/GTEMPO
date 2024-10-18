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



function split_mpotensor(mpoj::MPOTensor, trunc)
	ph = grassmannpspace()
	f = isomorphism(fuse(ph, ph), ph ⊗ ph)
	@tensor mpoj6[1,5,7;6,3,8] := mpoj[1,2,3,4] * conj(f[2,5,6]) * f[4,7,8]
	u, s, v = tsvd!(mpoj6, trunc=trunc)
	ss = sqrt(s)
	return permute(u * ss, (1,2), (4,3)), permute(ss * v, (1,2), (3,4))
end