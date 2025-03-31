function sysdynamics_imaginary!(gmps::FockMPS, lattice::AbstractFockLattice, model::AndersonIM; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a, b = siam_coeffs_fock(μ, U, -lattice.δτ)   
	for band in 1:lattice.bands
		for i in 1:lattice.Nτ
            pos = index(lattice, i, branch=:τ, band=band)
            apply!(exp(NTerm(pos, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end

	# interacting dynamics
	if U != zero(U)
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		for i in 1:lattice.Nτ
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i, branch=:τ, band=band)
				pos2 = index(lattice, i, branch=:τ, band=band+1)
				apply!(exp(NTerm(pos1, pos2, coeff=b)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end
	end
	return gmps
end
function sysdynamics_forward!(gmps::FockMPS, lattice::AbstractFockLattice, model::AndersonIM; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a, b = siam_coeffs_fock(μ, U, -im*lattice.δt) 
	for band in 1:lattice.bands
		for i in 1:lattice.Nt
            pos = index(lattice, i, branch=:+, band=band)
            apply!(exp(NTerm(pos, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end

	# interacting dynamics
	if U != zero(U)
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		for i in 1:lattice.Nt
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i, branch=:+, band=band)
				pos2 = index(lattice, i, branch=:+, band=band+1)
				apply!(exp(NTerm(pos1, pos2, coeff=b)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end
	end
	return gmps
end
function sysdynamics_backward!(gmps::FockMPS, lattice::AbstractFockLattice, model::AndersonIM; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a, b = siam_coeffs_fock(μ, U, im*lattice.δt) 
	for band in 1:lattice.bands
		for i in 1:lattice.Nt
			pos = index(lattice, i, branch=:-, band=band)
			apply!(exp(NTerm(pos, coeff=a)), gmps)
			canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	
		end
	end

	# interacting dynamics
	if U != 0.
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		for i in 1:lattice.Nt
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i, branch=:-, band=band)
				pos2 = index(lattice, i, branch=:-, band=band+1)
				apply!(exp(NTerm(pos1, pos2, coeff=b)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end	
	end
	return gmps
end

function siam_coeffs_fock(μ, U, dt)
	# a = exp(dt*(μ - U/2))
	a = exp(dt*μ) - 1
	b = exp(dt*U) - 1
	return a, b
end
