"""
	struct IRLM{L<:AbstractFermionicBath, R<:AbstractFermionicBath}

Interacting resonant level model, a spinless fermionic model with 
three impurities and two baths
"""
struct IRLM <: AbstractImpurityHamiltonian
	μ::Float64	
	J::Float64
	U::Float64
end
IRLM(; μ::Real, J::Real, U::Real) = IRLM(convert(Float64, μ), convert(Float64, J), convert(Float64, U))



# function sysdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::IRLM; trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	@assert lattice.bands == 3
# 	μ, J, U = model.μ, model.J, model.U
# 	δt = lattice.δt
# 	# forward evolution
# 	## ε part
# 	a = exp(-im*δt*(μ-U))
# 	for i in 1:lattice.k-1
# 		pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=2), index(lattice, i, conj=false, branch=:+, band=2)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=1), index(lattice, i, conj=false, branch=:+, band=1)
# 		apply!(exp(GTerm(pos1, pos2, coeff=1)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=3), index(lattice, i, conj=false, branch=:+, band=3)
# 		apply!(exp(GTerm(pos1, pos2, coeff=1)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
# 	end

# 	## J part
# 	a = -im*δt*J
# 	for i in 1:lattice.k-1
# 		pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=2), index(lattice, i, conj=false, branch=:+, band=1)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=1), index(lattice, i, conj=false, branch=:+, band=2)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=2), index(lattice, i, conj=false, branch=:+, band=3)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=3), index(lattice, i, conj=false, branch=:+, band=2)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
# 	end

# 	## U part
# 	a = -im*δt*U
# 	for i in 1:lattice.k-1
# 		pos1 = index(lattice, i+1, conj=true, branch=:+, band=1)
# 		pos2 = index(lattice, i+1, conj=true, branch=:+, band=2)
# 		pos3 = index(lattice, i, conj=false, branch=:+, band=2)
# 		pos4 = index(lattice, i, conj=false, branch=:+, band=1)
# 		apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	

# 		pos1 = index(lattice, i+1, conj=true, branch=:+, band=3)
# 		pos2 = index(lattice, i+1, conj=true, branch=:+, band=2)
# 		pos3 = index(lattice, i, conj=false, branch=:+, band=2)
# 		pos4 = index(lattice, i, conj=false, branch=:+, band=3)
# 		apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))					
# 	end

# 	# backward evolution
# 	## ε part
# 	a = exp(im*δt*(μ-U))
# 	for i in 1:lattice.k-1
# 		pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=2), index(lattice, i+1, conj=false, branch=:-, band=2)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=1), index(lattice, i+1, conj=false, branch=:-, band=1)
# 		apply!(exp(GTerm(pos1, pos2, coeff=1)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=3), index(lattice, i+1, conj=false, branch=:-, band=3)
# 		apply!(exp(GTerm(pos1, pos2, coeff=1)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
# 	end

# 	## J part
# 	a = im*δt*J
# 	for i in 1:lattice.k-1
# 		pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=2), index(lattice, i+1, conj=false, branch=:-, band=1)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=1), index(lattice, i+1, conj=false, branch=:-, band=2)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=2), index(lattice, i+1, conj=false, branch=:-, band=3)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

# 		pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=3), index(lattice, i+1, conj=false, branch=:-, band=2)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
# 	end

# 	## U part
# 	a = im*δt*U
# 	for i in 1:lattice.k-1
# 		pos1 = index(lattice, i, conj=true, branch=:-, band=1)
# 		pos2 = index(lattice, i, conj=true, branch=:-, band=2)
# 		pos3 = index(lattice, i+1, conj=false, branch=:-, band=2)
# 		pos4 = index(lattice, i+1, conj=false, branch=:-, band=1)
# 		apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	

# 		pos1 = index(lattice, i, conj=true, branch=:-, band=3)
# 		pos2 = index(lattice, i, conj=true, branch=:-, band=2)
# 		pos3 = index(lattice, i+1, conj=false, branch=:-, band=2)
# 		pos4 = index(lattice, i+1, conj=false, branch=:-, band=3)
# 		apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))					
# 	end

# 	return gmps
# end

sysdynamics_imaginary!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::IRLM; trunc::TruncationScheme=DMRG.DefaultTruncation) = sysdynamics_util!(
					gmps, lattice, model, lattice.Nτ, -lattice.δτ, :τ, trunc)
sysdynamics_forward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::IRLM; trunc::TruncationScheme=DMRG.DefaultTruncation) = sysdynamics_util!(
					gmps, lattice, model, lattice.Nt, -im*lattice.δt, :+, trunc)
sysdynamics_backward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::IRLM; trunc::TruncationScheme=DMRG.DefaultTruncation) = sysdynamics_util!(
					gmps, lattice, model, lattice.Nt, im*lattice.δt, :-, trunc)


function sysdynamics_util!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::IRLM, N::Int, dt::Number, branch::Symbol, trunc)
	@assert lattice.bands == 3
	μ, J, U = model.μ, model.J, model.U
	# forward evolution
	## ε part
	a = exp(dt*(μ-U))
	for i in 1:N
		x, y = (branch == :-) ? (i, i+1) : (i+1, i)
		pos1, pos2 = index(lattice, x, conj=true, branch=branch, band=2), index(lattice, y, conj=false, branch=branch, band=2)
		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

		pos1, pos2 = index(lattice, x, conj=true, branch=branch, band=1), index(lattice, y, conj=false, branch=branch, band=1)
		apply!(exp(GTerm(pos1, pos2, coeff=1)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

		pos1, pos2 = index(lattice, x, conj=true, branch=branch, band=3), index(lattice, y, conj=false, branch=branch, band=3)
		apply!(exp(GTerm(pos1, pos2, coeff=1)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
	end

	## J part
	a = dt*J
	for i in 1:N
		x, y = (branch == :-) ? (i, i+1) : (i+1, i)
		pos1, pos2 = index(lattice, x, conj=true, branch=branch, band=2), index(lattice, y, conj=false, branch=branch, band=1)
		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

		pos1, pos2 = index(lattice, x, conj=true, branch=branch, band=1), index(lattice, y, conj=false, branch=branch, band=2)
		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

		pos1, pos2 = index(lattice, x, conj=true, branch=branch, band=2), index(lattice, y, conj=false, branch=branch, band=3)
		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))

		pos1, pos2 = index(lattice, x, conj=true, branch=branch, band=3), index(lattice, y, conj=false, branch=branch, band=2)
		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
	end

	## U part
	a = dt*U
	for i in 1:N
		x, y = (branch == :-) ? (i, i+1) : (i+1, i)
		pos1 = index(lattice, x, conj=true, branch=branch, band=1)
		pos2 = index(lattice, x, conj=true, branch=branch, band=2)
		pos3 = index(lattice, y, conj=false, branch=branch, band=2)
		pos4 = index(lattice, y, conj=false, branch=branch, band=1)
		apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	

		pos1 = index(lattice, x, conj=true, branch=branch, band=3)
		pos2 = index(lattice, x, conj=true, branch=branch, band=2)
		pos3 = index(lattice, y, conj=false, branch=branch, band=2)
		pos4 = index(lattice, y, conj=false, branch=branch, band=3)
		apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))					
	end
	return gmps
end


# function hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, model::IRLM;
# 	leftcorr::Union{Nothing, <:RealCorrelationFunction}=nothing,
# 	rightcorr::Union{Nothing, <:RealCorrelationFunction}=nothing, kwargs...)
# 	if isnothing(leftcorr)
# 		leftcorr = correlationfunction(model.leftbath, lattice)
# 	end
# 	if isnothing(rightcorr)
# 		rightcorr = correlationfunction(model.rightbath, lattice)
# 	end
# 	gmps = hybriddynamics(gmps, lattice, leftcorr; band=1, kwargs...)
# 	gmps = hybriddynamics(gmps, lattice, leftcorr; band=3, kwargs...)
# 	return gmps
# end
