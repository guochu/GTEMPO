const DefaultKTruncation = truncdimcutoff(D=1000, ϵ=1.0e-10, add_back=0)

struct SISB{B <: AbstractFermionicBath} <: AbstractImpurityModel
	bath::B
	U::Float64
	μ::Float64
end
sys_size(x::SISB) = 1
SISB(bath::AbstractFermionicBath; U::Real, μ::Real) = SISB(bath, convert(Float64, U), convert(Float64, μ))

struct SIDB{L <: AbstractFermionicBath, R <: AbstractFermionicBath} <: AbstractImpurityModel
	leftbath::L
	rightbath::R
	U::Float64
	μ::Float64
end
sys_size(x::SIDB) = 1
SIDB(leftbath::AbstractFermionicBath, rightbath::AbstractFermionicBath; U::Real, μ::Real) = SIDB(
	leftbath, rightbath, convert(Float64, U), convert(Float64, μ))


const SingleImpurityModel = Union{SISB, SIDB}

function hybriddynamics(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, model::SISB; 
					corr::Union{Nothing, <:ImagCorrelationFunction}=nothing, trunc::TruncationScheme=DefaultITruncation)
	@assert lattice.β == model.bath.β
	if isnothing(corr)
		corr = correlationfunction(model.bath, lattice)
	end
	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
end

function hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SISB; 
					corr::Union{Nothing, <:RealCorrelationFunction}=nothing, 
					trunc::TruncationScheme=DefaultITruncation)
	if isnothing(corr)
		corr = correlationfunction(model.bath, lattice)
	end
	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
end


function hybriddynamics(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, model::SIDB; 
					corr::Union{Nothing, <:ImagCorrelationFunction}=nothing, trunc::TruncationScheme=DefaultITruncation)
	if isnothing(corr)
		corr = correlationfunction(model.leftbath, lattice) + correlationfunction(model.rightbath, lattice)
	end
	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
end

function hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, model::SIDB; 
					corr::Union{Nothing, <:RealCorrelationFunction}=nothing, 
					trunc::TruncationScheme=DefaultITruncation)
	if isnothing(corr)
		corr = correlationfunction(model.leftbath, lattice) + correlationfunction(model.rightbath, lattice)
	end
	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
end



# # # impurity dynamics
# """
# 	si_freedynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order; μ::Real, trunc::TruncationScheme=DMRG.DefaultTruncation)

# Imaginary time free dynamics of single impurity model
# """
# function si_freedynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order; μ::Real, trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	μ = convert(Float64, μ)
# 	a = exp(-lattice.δτ*μ)
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
#             pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
#             apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
#             canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 		end
# 	end		
# 	return gmps
# end
# si_freedynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...) = si_freedynamics!(copy(gmps), lattice; kwargs...)
# si_freedynamics(lattice::AbstractGrassmannLattice; kwargs...) = si_freedynamics(GrassmannMPS(eltype(lattice), length(lattice)), lattice; kwargs...)

# """
# 	si_interactiondynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order; U::Real, trunc::TruncationScheme=DMRG.DefaultTruncation)

# Imaginary time interaction dynamics of single impurity model
# """
# function si_interactiondynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order; U::Real, trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	U = convert(Float64, U)
# 	if U != 0.
# 		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
# 		a = -lattice.δτ*U
# 		for i in 1:lattice.k-1
# 			for band in 1:2:lattice.bands
# 				pos1 = index(lattice, i+1, conj=true, band=band)
# 				pos2 = index(lattice, i+1, conj=true, band=band+1)
# 				pos3 = index(lattice, i, conj=false, band=band+1)
# 				pos4 = index(lattice, i, conj=false, band=band)
# 				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
# 			end
# 		end
# 	end
# 	return gmps
# end
# si_interactiondynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...) = si_interactiondynamics!(copy(gmps), lattice; kwargs...)

# function si_sysdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice; 
# 						μ::Real, U::Real=0, trunc::TruncationScheme=DefaultKTruncation)
# 	# free dynamics
# 	μ = convert(Float64, μ)
# 	a = exp(-lattice.δτ*μ)
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
#             pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
#             apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
#             canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 		end
# 	end		

# 	# interacting dynamics
# 	U = convert(Float64, U)
# 	if U != 0.
# 		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
# 		# a = -lattice.δτ*U
# 		b = a^2 * (exp(-lattice.δτ*U) - 1)
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
# si_sysdynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...) = si_sysdynamics!(copy(gmps), lattice; kwargs...)
# si_sysdynamics(lattice::AbstractGrassmannLattice; kwargs...) = si_sysdynamics(GrassmannMPS(eltype(lattice), length(lattice)), lattice; kwargs...)
function sysdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation) 
	μ, U = model.μ, model.U
	a = exp(-lattice.δτ*μ)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
            pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end		

	# interacting dynamics
	U = convert(Float64, U)
	if U != 0.
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -lattice.δτ*U
		b = a^2 * (exp(-lattice.δτ*U) - 1)
		for i in 1:lattice.k-1
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i+1, conj=true, band=band)
				pos2 = index(lattice, i+1, conj=true, band=band+1)
				pos3 = index(lattice, i, conj=false, band=band+1)
				pos4 = index(lattice, i, conj=false, band=band)
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
			end
		end
	end
	return gmps	
end

# """
# 	si_sysdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order; μ::Real, trunc::TruncationScheme=DMRG.DefaultTruncation)

# Real time free dynamics of single impurity model
# """
# function si_freedynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice; 
# 						μ::Real, sys_state::Symbol=:vacuum, trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	(sys_state in (:filled, :vacuum)) || throw(ArgumentError("sys_state must be :filled, :vacuum, or :thermal"))
# 	a = exp(-im*lattice.δt*μ)
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
#             pos1, pos2 = index(lattice, i+1, conj=true, forward=true, band=band), index(lattice, i, conj=false, forward=true, band=band)
#             apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
#             canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 		end
# 	end
# 	# if sys_state == :thermal
# 	# 	a = exp(lattice.β*μ)
# 	# 	for band in 1:lattice.bands
# 	# 		pos1, pos2 = index(lattice, 1, conj=true, forward=true, band=band), index(lattice, 1, conj=false, forward=false, band=band)
# 	# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 	# 		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	
# 	# 	end
# 	# elseif sys_state == :filled
# 	# 	error("not implemented for sys_state=:filled")
# 	# end
# 	(sys_state == :filled) && throw(ArgumentError("not implemented for sys_state=:filled"))
# 	a = exp(im*lattice.δt*μ)
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
# 			pos1, pos2 = index(lattice, i, conj=true, forward=false, band=band), index(lattice, i+1, conj=false, forward=false, band=band)
# 			apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 			canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	
# 		end
# 	end
# 	return gmps
# end

# function si_interactiondynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice; U::Real, trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	U = convert(Float64, U)
# 	if U != 0.
# 		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
# 		a = -im*lattice.δt*U
# 		for i in 1:lattice.k-1
# 			for band in 1:2:lattice.bands
# 				pos1 = index(lattice, i+1, conj=true, forward=true, band=band)
# 				pos2 = index(lattice, i+1, conj=true, forward=true, band=band+1)
# 				pos3 = index(lattice, i, conj=false, forward=true, band=band+1)
# 				pos4 = index(lattice, i, conj=false, forward=true, band=band)
# 				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 			end
# 		end
# 		a = im*lattice.δt*U
# 		for i in 1:lattice.k-1
# 			for band in 1:2:lattice.bands
# 				pos1 = index(lattice, i, conj=true, forward=false, band=band)
# 				pos2 = index(lattice, i, conj=true, forward=false, band=band+1)
# 				pos3 = index(lattice, i+1, conj=false, forward=false, band=band+1)
# 				pos4 = index(lattice, i+1, conj=false, forward=false, band=band)
# 				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
# 			end
# 		end	
# 	end
# 	return gmps
# end

function sysdynamics_forward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a = exp(-im*lattice.δt*μ)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
            pos1, pos2 = index(lattice, i+1, conj=true, forward=true, band=band), index(lattice, i, conj=false, forward=true, band=band)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end

	# interacting dynamics
	if U != 0.
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		b = a^2 * (exp(-im*lattice.δt*U) - 1)
		for i in 1:lattice.k-1
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i+1, conj=true, forward=true, band=band)
				pos2 = index(lattice, i+1, conj=true, forward=true, band=band+1)
				pos3 = index(lattice, i, conj=false, forward=true, band=band+1)
				pos4 = index(lattice, i, conj=false, forward=true, band=band)
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end
	end
	return gmps
end
function sysdynamics_backward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a = exp(-im*lattice.δt*μ)
	ac = conj(a)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
			pos1, pos2 = index(lattice, i, conj=true, forward=false, band=band), index(lattice, i+1, conj=false, forward=false, band=band)
			apply!(exp(GTerm(pos1, pos2, coeff=ac)), gmps)
			canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))	
		end
	end

	# interacting dynamics
	if U != 0.
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		b = a^2 * (exp(-im*lattice.δt*U) - 1)
		bc = conj(b)
		for i in 1:lattice.k-1
			for band in 1:2:lattice.bands
				pos1 = index(lattice, i, conj=true, forward=false, band=band)
				pos2 = index(lattice, i, conj=true, forward=false, band=band+1)
				pos3 = index(lattice, i+1, conj=false, forward=false, band=band+1)
				pos4 = index(lattice, i+1, conj=false, forward=false, band=band)
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=bc)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end	
	end
	return gmps
end

# function systhermalstate!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
# 	μ, U = model.μ, model.U
# 	β = model.bath.β
# 	a = exp(-β*μ)
# 	for band in 1:lattice.bands
# 		pos1, pos2 = index(lattice, 1, conj=true, forward=true, band=band), index(lattice, 1, conj=false, forward=false, band=band)
# 		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		canonicalize!(gmps, alg=Orthogonalize(trunc=trunc))	
# 	end	
# 	if U != 0.
# 		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
# 		b = a^2 * (exp(-β*U) - 1)
# 		for band in 1:2:lattice.bands
# 			pos1 = index(lattice, 1, conj=true, forward=true, band=band)
# 			pos2 = index(lattice, 1, conj=true, forward=true, band=band+1)
# 			pos3 = index(lattice, 1, conj=false, forward=false, band=band+1)
# 			pos4 = index(lattice, 1, conj=false, forward=false, band=band)
# 			apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
# 			canonicalize!(gmps, alg=Orthogonalize(trunc=trunc))			
# 		end
# 	end
# 	return gmps
# end

function si_sysdynamics_stepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice; μ::Real, U::Real=0, trunc::Union{Nothing,TruncationScheme}=DefaultKTruncation)
	# free dynamics
	a = exp(-im*lattice.δt*μ)
	i = lattice.k - 1
	for band in 1:lattice.bands
		pos1, pos2 = index(lattice, i+1, conj=true, forward=true, band=band), index(lattice, i, conj=false, forward=true, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
	end
	ac = conj(a)
	for band in 1:lattice.bands
		pos1, pos2 = index(lattice, i, conj=true, forward=false, band=band), index(lattice, i+1, conj=false, forward=false, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=ac)), gmps)
	end	

	# interaction dynamics
	U = convert(Float64, U)
	if U != 0.
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		b = a^2 * (exp(-im*lattice.δt*U) - 1)
		for band in 1:2:lattice.bands
			pos1 = index(lattice, i+1, conj=true, forward=true, band=band)
			pos2 = index(lattice, i+1, conj=true, forward=true, band=band+1)
			pos3 = index(lattice, i, conj=false, forward=true, band=band+1)
			pos4 = index(lattice, i, conj=false, forward=true, band=band)
			apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
		end
		# a = im*lattice.δt*U
		bc = conj(b)
		for band in 1:2:lattice.bands
			pos1 = index(lattice, i, conj=true, forward=false, band=band)
			pos2 = index(lattice, i, conj=true, forward=false, band=band+1)
			pos3 = index(lattice, i+1, conj=false, forward=false, band=band+1)
			pos4 = index(lattice, i+1, conj=false, forward=false, band=band)
			apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=bc)), gmps)
		end
	end	
	canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
	return gmps
end
sysdynamicsstepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; kwargs...) = si_sysdynamics_stepper!(gmps, lattice; μ=model.μ, U=model.U, kwargs...)

# function si_sysdynamics(lattice::RealGrassmannLattice; μ::Real, U::Real=0, sys_state::Symbol=:vacuum, trunc::TruncationScheme=DMRG.DefaultTruncation)
# 	# free dynamics
# 	gmps = vacuumstate(lattice)
# 	(sys_state in (:filled, :vacuum)) || throw(ArgumentError("sys_state must be :filled, :vacuum, or :thermal"))
# 	a = exp(-im*lattice.δt*μ)
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
#             pos1, pos2 = index(lattice, i+1, conj=true, forward=true, band=band), index(lattice, i, conj=false, forward=true, band=band)
#             apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		end
# 	end
# 	(sys_state == :filled) && throw(ArgumentError("not implemented for sys_state=:filled"))
# 	a = exp(im*lattice.δt*μ)
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
# 			pos1, pos2 = index(lattice, i, conj=true, forward=false, band=band), index(lattice, i+1, conj=false, forward=false, band=band)
# 			apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 		end
# 	end

# 	# interaction dynamics
# 	U = convert(Float64, U)
# 	if U != 0.
# 		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
# 		a = -im*lattice.δt*U
# 		for i in 1:lattice.k-1
# 			for band in 1:2:lattice.bands
# 				pos1 = index(lattice, i+1, conj=true, forward=true, band=band)
# 				pos2 = index(lattice, i+1, conj=true, forward=true, band=band+1)
# 				pos3 = index(lattice, i, conj=false, forward=true, band=band+1)
# 				pos4 = index(lattice, i, conj=false, forward=true, band=band)
# 				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 			end
# 		end
# 		a = im*lattice.δt*U
# 		for i in 1:lattice.k-1
# 			for band in 1:2:lattice.bands
# 				pos1 = index(lattice, i, conj=true, forward=false, band=band)
# 				pos2 = index(lattice, i, conj=true, forward=false, band=band+1)
# 				pos3 = index(lattice, i+1, conj=false, forward=false, band=band+1)
# 				pos4 = index(lattice, i+1, conj=false, forward=false, band=band)
# 				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
# 			end
# 		end	
# 	end	
# 	return gmps
# end
# sysdynamics(lattice::AbstractGrassmannLattice, model::SingleImpurityModel; kwargs...) = si_sysdynamics(lattice; μ=model.μ, U=model.U, kwargs...)


