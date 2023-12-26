

function _leftorth!(psi::GrassmannMPS, alg::QR, trunc::TruncationScheme)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in 1:L-1
		q, r = leftorth!(psi[i], alg = alg)
		psi[i] = q
		nr = norm(r)
		_rescaling!(psi, nr)
		r = rmul!(r, 1/nr)
		psi[i + 1] = @tensor tmp[1 3; 4] := r[1,2] * psi[i+1][2,3,4]
	end
	return psi
end

function _rightorth!(psi::GrassmannMPS, alg::QR, trunc::TruncationScheme)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in L:-1:2
		l, q = rightorth(psi[i], (1,), (2, 3), alg=LQ())
		psi[i] = permute(q, (1,2), (3,))
		nl = norm(l)
		_rescaling!(psi, nl)
		l = rmul!(l, 1/nl)
		psi[i-1] = @tensor tmp[1 2; 4] := psi[i-1][1,2,3] * l[3,4] 
	end
	return psi
end

function _rightorth!(psi::GrassmannMPS, alg::SVD, trunc::TruncationScheme)
	L = length(psi)
	for i in L:-1:2
		u, s, v, err = DMRG.stable_tsvd(psi[i], (1,), (2, 3), trunc=trunc)
		psi[i] = permute(v, (1,2), (3,))
		u2 = u * s
		nl = norm(u2)
		_rescaling!(psi, nl)
		u2 = rmul!(u2, 1/nl)
		psi[i-1] = @tensor tmp[-1 -2; -3] := psi[i-1][-1, -2, 1] * u2[1, -3]
	end
	return psi
end

function DMRG.canonicalize!(psi::GrassmannMPS; alg::Orthogonalize = Orthogonalize(SVD(), DMRG.DefaultTruncation, normalize=false))
	alg.normalize && @warn "canonicalize with normalization not implemented for GrassmannMPS"
	L = length(psi)
	_leftorth!(psi, QR(), NoTruncation())
	_rightorth!(psi, alg.orth, alg.trunc)

	return _rescaling!(psi)
end

function _rescaling!(psi::GrassmannMPS, n::Real)
	L = length(psi)
	scale1 = n^(1/L)
	psi.scale[] = scale(psi) * scale1
	return psi
end
function _rescaling!(psi::GrassmannMPS)
	nrm1 = norm(psi[1])
	psi[1] = rmul!(psi[1], 1/nrm1)
	return _rescaling!(psi, nrm1)
end
# function _rescaling!(psi::GrassmannMPS)
# 	nrm1 = norm(psi[1])
# 	L = length(psi)
# 	scale1 = nrm1^(1/L)
# 	psi[1] = rmul!(psi[1], 1/nrm1)
# 	psi.scale[] = scale(psi) * scale1
# 	return psi
# end