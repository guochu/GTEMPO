
TK.leftorth!(psi::GrassmannMPS; alg::Orthogonalize = Orthogonalize()) = _leftorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _leftorth!(psi::GrassmannMPS, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in 1:L-1
		q, r = leftorth!(GrassmannTensorMap(psi[i]), alg = alg)
		psi[i] = get_data(q)
		# nr = norm(r)
		# (nr ≈ zero(nr)) && @warn "norm of GrassmannMPS is zero"
		# _rescaling!(psi, nr)
		# r = rmul!(r, 1/nr)
		_renormalize!(psi, get_data(r), normalize)
		@tensor tmp[1 3; 4] := r[1,2] * GrassmannTensorMap(psi[i+1])[2,3,4]
		psi[i + 1] = get_data(tmp)
	end
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function _leftorth!(psi::GrassmannMPS, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	L = length(psi)
	# errs = Float64[]
	maxerr = 0.
	for i in 1:L-1
		u, s, v, err = stable_tsvd!(GrassmannTensorMap(psi[i]), trunc=trunc)
		nr = _renormalize!(psi, get_data(s), normalize)
		rerror = sqrt(err * err / (nr * nr + err * err))
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", rerror)
		psi[i] = get_data(u)
		v2 = s * v
		@tensor tmp[-1 -2; -3] := v2[-1, 1] * GrassmannTensorMap(psi[i+1])[1,-2,-3]
		psi[i+1] = get_data(tmp)
		psi.s[i+1] = get_data(s)
		# push!(errs, err)
		maxerr = max(maxerr, rerror)
	end
	(verbosity > 0) && println("Max SVD truncerror in leftorth: ", maxerr)
	_renormalize!(psi, psi[L], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

TK.rightorth!(psi::GrassmannMPS; alg::Orthogonalize = Orthogonalize()) = _rightorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
function _rightorth!(psi::GrassmannMPS, alg::QR, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	!isa(trunc, NoTruncation) &&  @warn "truncation has no effect with QR"
	L = length(psi)
	for i in L:-1:2
		l, q = rightorth(GrassmannTensorMap(psi[i]), (1,), (2, 3), alg=LQ())
		psi[i] = get_data(permute(q, (1,2), (3,)))
		# nl = norm(l)
		# (nl ≈ zero(nl)) && @warn "norm of GrassmannMPS is zero"
		# _rescaling!(psi, nl)
		# l = rmul!(l, 1/nl)
		_renormalize!(psi, get_data(l), normalize)
		@tensor tmp[1 2; 4] := GrassmannTensorMap(psi[i-1])[1,2,3] * l[3,4] 
		psi[i-1] = get_data(tmp)
	end
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function _rightorth!(psi::GrassmannMPS, alg::SVD, trunc::TruncationScheme, normalize::Bool, verbosity::Int)
	L = length(psi)
	maxerr = 0.
	for i in L:-1:2
		u, s, v, err = stable_tsvd(GrassmannTensorMap(psi[i]), (1,), (2, 3), trunc=trunc)
		psi[i] = get_data(permute(v, (1,2), (3,)))
		nr = _renormalize!(psi, get_data(s), normalize)
		rerror = sqrt(err * err / (nr * nr + err * err))
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", rerror)
		u2 = u * s
		# nl = norm(u2)
		# (nl ≈ zero(nl)) && @warn "norm of GrassmannMPS is zero"
		# _rescaling!(psi, nl)
		# u2 = rmul!(u2, 1/nl)
		@tensor tmp[-1 -2; -3] := GrassmannTensorMap(psi[i-1])[-1, -2, 1] * u2[1, -3]
		psi[i-1] = get_data(tmp)
		psi.s[i] = get_data(s)
		maxerr = max(maxerr, err)
	end
	(verbosity > 0) && println("Max SVD truncerror in rightorth: ", maxerr)
	_renormalize!(psi, psi[1], normalize)
	_renormalize_coeff!(psi, normalize)
	return psi
end

function DMRG.canonicalize!(psi::GrassmannMPS; alg::Orthogonalize = Orthogonalize(trunc=DMRG.DefaultTruncation, normalize=false))
	alg.normalize && @warn "canonicalize with renormalization not recommanded for GrassmannMPS"
	L = length(psi)
	_leftorth!(psi, QR(), NoTruncation(), alg.normalize, alg.verbosity)
	_rightorth!(psi, alg.orth, alg.trunc, alg.normalize, alg.verbosity)
	return psi
end

function _rescaling!(psi, n::Real)
	L = length(psi)
	scale1 = n^(1/L)
	setscaling!(psi, scaling(psi) * scale1)
	return psi
end
function _rescaling!(psi)
	nrm1 = norm(psi[1])
	psi[1] = rmul!(psi[1], 1/nrm1)
	return _rescaling!(psi, nrm1)
end


function _renormalize!(psi, r, normalize::Bool)
	nr = norm(r)
	if nr != zero(nr)
		if !normalize
			_rescaling!(psi, nr)
		end
		lmul!(1/nr, r)  
  	end
  	return nr
end

function _renormalize_coeff!(psi, normalize::Bool)
	normalize && setscaling!(psi, 1)
end
