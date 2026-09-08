"""
	struct AndersonIM{B <: AbstractFermionicBath}

Single-orbital Anderson impurity model with one bath
"""
struct AndersonIM <: AbstractImpurityHamiltonian
	U::Float64
	μ::Float64
end
AndersonIM(; U::Real, μ::Real) = AndersonIM(convert(Float64, U), convert(Float64, μ))

# function hybriddynamics(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, model::AndersonIM; 
# 					corr::Union{Nothing, <:ImagCorrelationFunction}=nothing, trunc::TruncationScheme=DefaultITruncation)
# 	@assert lattice.β == model.bath.β
# 	if isnothing(corr)
# 		corr = correlationfunction(model.bath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end

# function hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AndersonIM; 
# 					corr::Union{Nothing, <:RealCorrelationFunction}=nothing, 
# 					trunc::TruncationScheme=DefaultITruncation)
# 	if isnothing(corr)
# 		corr = correlationfunction(model.bath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end

# function hybriddynamics(gmps::GrassmannMPS, lattice::MixedGrassmannLattice, model::AndersonIM; 
# 					corr::Union{Nothing, <:MixedCorrelationFunction}=nothing, 
# 					trunc::TruncationScheme=DefaultITruncation)
# 	@assert lattice.β == model.bath.β
# 	if isnothing(corr)
# 		corr = correlationfunction(model.bath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end

# function hybriddynamics(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, model::SIDB; 
# 					corr::Union{Nothing, <:ImagCorrelationFunction}=nothing, trunc::TruncationScheme=DefaultITruncation)
# 	if isnothing(corr)
# 		corr = correlationfunction(model.leftbath, lattice) + correlationfunction(model.rightbath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end

# function hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, model::SIDB; 
# 					corr::Union{Nothing, <:RealCorrelationFunction}=nothing, 
# 					trunc::TruncationScheme=DefaultITruncation)
# 	if isnothing(corr)
# 		corr = correlationfunction(model.leftbath, lattice) + correlationfunction(model.rightbath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end


# function sysdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation) 
# 	μ, U = model.μ, model.U
# 	a, b = siam_coeffs(μ, U, -lattice.δτ)
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
#             pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
#             apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
#             canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 		end
# 	end		

# 	# interacting dynamics
# 	if U != zero(U)
# 		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
# 		# a = -lattice.δτ*U
# 		# b = a^2 * (exp(-lattice.δτ*U) - 1)
# 		for i in 1:lattice.k-1
# 			for band in 1:2:lattice.bands
# 				pos1 = index(lattice, i+1, conj=true, band=band)
# 				pos2 = index(lattice, i+1, conj=true, band=band+1)
# 				pos3 = index(lattice, i, conj=false, band=band+1)
# 				pos4 = index(lattice, i, conj=false, band=band)
# 				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
# 				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
# 			end
# 		end
# 	end
# 	return gmps	
# end


# function sysdynamics_forward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
# 	# free dynamics
# 	μ, U = model.μ, model.U
# 	a, b = siam_coeffs(μ, U, -im*lattice.δt) 
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
#             pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=band), index(lattice, i, conj=false, branch=:+, band=band)
#             apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
#             canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 		end
# 	end

# 	# interacting dynamics
# 	if U != zero(U)
# 		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
# 		# a = -im*lattice.δt*U
# 		# b = a^2 * (exp(-im*lattice.δt*U) - 1)
# 		for i in 1:lattice.k-1
# 			for band in 1:2:lattice.bands
# 				pos1 = index(lattice, i+1, conj=true, branch=:+, band=band)
# 				pos2 = index(lattice, i+1, conj=true, branch=:+, band=band+1)
# 				pos3 = index(lattice, i, conj=false, branch=:+, band=band+1)
# 				pos4 = index(lattice, i, conj=false, branch=:+, band=band)
# 				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
# 				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 			end
# 		end
# 	end
# 	return gmps
# end
# function sysdynamics_backward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
# 	# free dynamics
# 	μ, U = model.μ, model.U
# 	# a = exp(-im*lattice.δt*μ)
# 	# ac = conj(a)
# 	a, b = siam_coeffs(μ, U, im*lattice.δt) 
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
# 			pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=band), index(lattice, i+1, conj=false, branch=:-, band=band)
# 			apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 			canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	
# 		end
# 	end

# 	# interacting dynamics
# 	if U != 0.
# 		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
# 		# a = -im*lattice.δt*U
# 		# b = a^2 * (exp(-im*lattice.δt*U) - 1)
# 		# bc = conj(b)
# 		for i in 1:lattice.k-1
# 			for band in 1:2:lattice.bands
# 				pos1 = index(lattice, i, conj=true, branch=:-, band=band)
# 				pos2 = index(lattice, i, conj=true, branch=:-, band=band+1)
# 				pos3 = index(lattice, i+1, conj=false, branch=:-, band=band+1)
# 				pos4 = index(lattice, i+1, conj=false, branch=:-, band=band)
# 				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
# 				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 			end
# 		end	
# 	end
# 	return gmps
# end

"""
	sysdynamics_imaginary!(gmps, lattice, model::AndersonIM; trunc) -> GrassmannMPS

Analytical solution of the imaginary-time propagator for the Anderson
impurity: the Hamiltonian terms μ(n₁+n₂) and U n₁n₂ all commute, so the
factorization into GTerms is exact. Building an equivalent model from
`ImpurityHamiltonian` and calling the generic `sysdynamics_imaginary!`
yields exactly the same result.
"""
function sysdynamics_imaginary!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AndersonIM; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a, b = siam_coeffs(μ, U, -lattice.δτ)   
	for band in 1:lattice.bands
		for i in 1:lattice.kτ-1
            pos1, pos2 = index(lattice, i+1, conj=true, branch=:τ, band=band), index(lattice, i, conj=false, branch=:τ, band=band)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end

	# interacting dynamics
	if U != zero(U)
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		# b = a^2 * (exp(-lattice.δτ*U) - 1)
		for i in 1:lattice.Nτ
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i+1, conj=true, branch=:τ, band=band)
				pos2 = index(lattice, i+1, conj=true, branch=:τ, band=band+1)
				pos3 = index(lattice, i, conj=false, branch=:τ, band=band+1)
				pos4 = index(lattice, i, conj=false, branch=:τ, band=band)
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end
	end
	return gmps
end
"""
	sysdynamics_forward!(gmps, lattice, model::AndersonIM; trunc) -> GrassmannMPS

Analytical solution of the forward-branch propagator for the Anderson
impurity: the Hamiltonian terms μ(n₁+n₂) and U n₁n₂ all commute, so the
factorization into GTerms is exact. Building an equivalent model from
`ImpurityHamiltonian` and calling the generic `sysdynamics_forward!`
yields exactly the same result.
"""
function sysdynamics_forward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AndersonIM; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a, b = siam_coeffs(μ, U, -im*lattice.δt) 
	for band in 1:lattice.bands
		for i in 1:lattice.Nt
            pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=band), index(lattice, i, conj=false, branch=:+, band=band)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end

	# interacting dynamics
	if U != zero(U)
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		# b = a^2 * (exp(-im*lattice.δt*U) - 1)
		for i in 1:lattice.Nt
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i+1, conj=true, branch=:+, band=band)
				pos2 = index(lattice, i+1, conj=true, branch=:+, band=band+1)
				pos3 = index(lattice, i, conj=false, branch=:+, band=band+1)
				pos4 = index(lattice, i, conj=false, branch=:+, band=band)
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end
	end
	return gmps
end
"""
	sysdynamics_backward!(gmps, lattice, model::AndersonIM; trunc) -> GrassmannMPS

Analytical solution of the backward-branch propagator for the Anderson
impurity: the Hamiltonian terms μ(n₁+n₂) and U n₁n₂ all commute, so the
factorization into GTerms is exact. Building an equivalent model from
`ImpurityHamiltonian` and calling the generic `sysdynamics_backward!`
yields exactly the same result.
"""
function sysdynamics_backward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model::AndersonIM; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	# a = exp(-im*lattice.δt*μ)
	# ac = conj(a)
	a, b = siam_coeffs(μ, U, im*lattice.δt) 
	for band in 1:lattice.bands
		for i in 1:lattice.Nt
			pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=band), index(lattice, i+1, conj=false, branch=:-, band=band)
			apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
			canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	
		end
	end

	# interacting dynamics
	if U != 0.
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		# b = a^2 * (exp(-im*lattice.δt*U) - 1)
		# bc = conj(b)
		for i in 1:lattice.Nt
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i, conj=true, branch=:-, band=band)
				pos2 = index(lattice, i, conj=true, branch=:-, band=band+1)
				pos3 = index(lattice, i+1, conj=false, branch=:-, band=band+1)
				pos4 = index(lattice, i+1, conj=false, branch=:-, band=band)
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end	
	end
	return gmps
end

function si_sysdynamics_stepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice; μ::Real, U::Real=0, trunc::Union{Nothing,TruncationScheme}=DefaultKTruncation)
	# free dynamics
	a = exp(-im*lattice.δt*μ)
	i = lattice.k - 1
	for band in 1:lattice.bands
		pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=band), index(lattice, i, conj=false, branch=:+, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
	end
	ac = conj(a)
	for band in 1:lattice.bands
		pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=band), index(lattice, i+1, conj=false, branch=:-, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=ac)), gmps)
	end	

	# interaction dynamics
	U = convert(Float64, U)
	if U != 0.
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		b = a^2 * (exp(-im*lattice.δt*U) - 1)
		for band in 1:2:lattice.bands
			pos1 = index(lattice, i+1, conj=true, branch=:+, band=band)
			pos2 = index(lattice, i+1, conj=true, branch=:+, band=band+1)
			pos3 = index(lattice, i, conj=false, branch=:+, band=band+1)
			pos4 = index(lattice, i, conj=false, branch=:+, band=band)
			apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
		end
		# a = im*lattice.δt*U
		bc = conj(b)
		for band in 1:2:lattice.bands
			pos1 = index(lattice, i, conj=true, branch=:-, band=band)
			pos2 = index(lattice, i, conj=true, branch=:-, band=band+1)
			pos3 = index(lattice, i+1, conj=false, branch=:-, band=band+1)
			pos4 = index(lattice, i+1, conj=false, branch=:-, band=band)
			apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=bc)), gmps)
		end
	end	
	canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
	return gmps
end
"""
    sysdynamicsstepper!(gmps, lattice::RealGrassmannLattice, model::AndersonIM; trunc=DefaultKTruncation)

Single-step version of `sysdynamics` for the Anderson impurity model: apply the
free (and, for two bands, the interaction) propagator of the current time step
`lattice.k - 1` onto `gmps` in place. Intended for online/stepwise evolutions
combined with `makestep` and `hybriddynamicsstepper!`.
"""
sysdynamicsstepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AndersonIM; kwargs...) = si_sysdynamics_stepper!(gmps, lattice; μ=model.μ, U=model.U, kwargs...)


function siam_coeffs(μ, U, dt)
	# a = exp(dt*(μ - U/2))
	a = exp(dt*μ)
	b = a^2 * (exp(dt*U) - 1)
	return a, b
end

"""
	systhermalstate!(gmps, lattice, model::AndersonIM; β, trunc) -> GrassmannMPS

Analytical solution of the normalized thermal state exp(-βĤ)/tr(exp(-βĤ))
for the Anderson impurity: the Hamiltonian terms μ(n₁+n₂) and U n₁n₂ all
commute, so the factorization of exp(-βĤ) into GTerms is exact. Building
an equivalent model from `ImpurityHamiltonian` and calling the generic
`systhermalstate!` yields exactly the same result. For `β == Inf` the
ground state projector is built via the generic Fock-space construction.
"""
function systhermalstate!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AndersonIM; β::Real, trunc::TruncationScheme=DefaultKTruncation)
	μ, U = model.μ, model.U
	if β == Inf
		return sysinitialstate!(gmps, lattice, fock_thermalstate(model, β, lattice.bands); trunc=trunc)
	end
	a, b = siam_coeffs(μ, U, -β)
	for band in 1:lattice.bands
		pos1, pos2 = index(lattice, 1, conj=true, branch=:+, band=band), index(lattice, 1, conj=false, branch=:-, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
	end
	if U != zero(U)
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		pos1 = index(lattice, 1, conj=true, branch=:+, band=1)
		pos2 = index(lattice, 1, conj=true, branch=:+, band=2)
		pos3 = index(lattice, 1, conj=false, branch=:-, band=2)
		pos4 = index(lattice, 1, conj=false, branch=:-, band=1)
		apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)			
	end
	# normalize: divide by tr(exp(-βĤ)) so that the result equals the
	# normalized thermal state, consistent with fock_thermalstate
	Z = (1 + a)^lattice.bands + (U != zero(U) ? b : 0)
	gmps[1] = gmps[1] * (1 / Z)
	canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
	return gmps
end 
