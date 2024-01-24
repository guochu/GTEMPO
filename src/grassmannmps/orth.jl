
TK.leftorth!(psi::GrassmannMPS; alg::Orthogonalize = Orthogonalize()) = _leftorth!(psi, alg.orth, alg.trunc, alg.normalize)
function _leftorth!(psi::GrassmannMPS, alg::QR, trunc::TruncationScheme, normalize::Bool)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in 1:L-1
		q, r = leftorth!(psi[i], alg = alg)
		psi[i] = q
		# nr = norm(r)
		# (nr ≈ zero(nr)) && @warn "norm of GrassmannMPS is zero"
		# _rescaling!(psi, nr)
		# r = rmul!(r, 1/nr)
		_renormalize!(psi, r, normalize)
		psi[i + 1] = @tensor tmp[1 3; 4] := r[1,2] * psi[i+1][2,3,4]
	end
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function _leftorth!(psi::GrassmannMPS, alg::SVD, trunc::TruncationScheme, normalize::Bool)
	L = length(psi)
	# errs = Float64[]
	for i in 1:L-1
		u, s, v, err = stable_tsvd!(psi[i], trunc=trunc)
		_renormalize!(psi, s, normalize)
		psi[i] = u
		v2 = s * v
		psi[i+1] = @tensor tmp[-1 -2; -3] := v2[-1, 1] * psi[i+1][1,-2,-3]
		psi.s[i+1] = s
		# push!(errs, err)
	end
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

TK.rightorth!(psi::GrassmannMPS; alg::Orthogonalize = Orthogonalize()) = _rightorth!(psi, alg.orth, alg.trunc, alg.normalize)
function _rightorth!(psi::GrassmannMPS, alg::QR, trunc::TruncationScheme, normalize::Bool)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in L:-1:2
		l, q = rightorth(psi[i], (1,), (2, 3), alg=LQ())
		psi[i] = permute(q, (1,2), (3,))
		# nl = norm(l)
		# (nl ≈ zero(nl)) && @warn "norm of GrassmannMPS is zero"
		# _rescaling!(psi, nl)
		# l = rmul!(l, 1/nl)
		_renormalize!(psi, l, normalize)
		psi[i-1] = @tensor tmp[1 2; 4] := psi[i-1][1,2,3] * l[3,4] 
	end
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function _rightorth!(psi::GrassmannMPS, alg::SVD, trunc::TruncationScheme, normalize::Bool)
	L = length(psi)
	for i in L:-1:2
		u, s, v, err = stable_tsvd(psi[i], (1,), (2, 3), trunc=trunc)
		psi[i] = permute(v, (1,2), (3,))
		_renormalize!(psi, s, normalize)
		u2 = u * s
		# nl = norm(u2)
		# (nl ≈ zero(nl)) && @warn "norm of GrassmannMPS is zero"
		# _rescaling!(psi, nl)
		# u2 = rmul!(u2, 1/nl)
		psi[i-1] = @tensor tmp[-1 -2; -3] := psi[i-1][-1, -2, 1] * u2[1, -3]
		psi.s[i] = s
	end
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function DMRG.canonicalize!(psi::GrassmannMPS; alg::Orthogonalize = Orthogonalize(trunc=DMRG.DefaultTruncation, normalize=false))
	alg.normalize && @warn "canonicalize with normalization not recommanded for GrassmannMPS"
	L = length(psi)
	_leftorth!(psi, QR(), NoTruncation(), alg.normalize)
	_rightorth!(psi, alg.orth, alg.trunc, alg.normalize)
	return psi
end

function _rescaling!(psi::GrassmannMPS, n::Real)
	L = length(psi)
	scale1 = n^(1/L)
	setscaling!(psi, scaling(psi) * scale1)
	return psi
end
function _rescaling!(psi::GrassmannMPS)
	nrm1 = norm(psi[1])
	psi[1] = rmul!(psi[1], 1/nrm1)
	return _rescaling!(psi, nrm1)
end


function _renormalize!(psi::GrassmannMPS, r, normalize::Bool)
	nr = norm(r)
	if nr != zero(nr)
		if !normalize
			_rescaling!(psi, nr)
		end
		lmul!(1/nr, r)  
  	end
end

function _renormalize_coeff!(psi::GrassmannMPS, normalize::Bool)
	normalize && setscaling!(psi, 1)
end
