"""
	struct SISB{B <: AbstractFermionicBath}

Single-orbital Anderson impurity model with one bath
"""
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
	if U != zero(U)
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


function sysdynamics_forward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a = exp(-im*lattice.δt*μ)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
            pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=band), index(lattice, i, conj=false, branch=:+, band=band)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end

	# interacting dynamics
	if U != zero(U)
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		b = a^2 * (exp(-im*lattice.δt*U) - 1)
		for i in 1:lattice.k-1
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
function sysdynamics_backward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a = exp(-im*lattice.δt*μ)
	ac = conj(a)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
			pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=band), index(lattice, i+1, conj=false, branch=:-, band=band)
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
				pos1 = index(lattice, i, conj=true, branch=:-, band=band)
				pos2 = index(lattice, i, conj=true, branch=:-, band=band+1)
				pos3 = index(lattice, i+1, conj=false, branch=:-, band=band+1)
				pos4 = index(lattice, i+1, conj=false, branch=:-, band=band)
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=bc)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
			end
		end	
	end
	return gmps
end

function sysdynamics_imaginary!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a = exp(-lattice.δτ*μ)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
            pos1, pos2 = index(lattice, i+1, conj=true, branch=:τ, band=band), index(lattice, i, conj=false, branch=:τ, band=band)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end

	# interacting dynamics
	if U != zero(U)
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		b = a^2 * (exp(-lattice.δτ*U) - 1)
		for i in 1:lattice.k-1
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
function sysdynamics_forward!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a = exp(-im*lattice.δt*μ)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
            pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=band), index(lattice, i, conj=false, branch=:+, band=band)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))			
		end
	end

	# interacting dynamics
	if U != zero(U)
		(lattice.bands == 2) || throw(ArgumentError("lattice should have two bands"))
		# a = -im*lattice.δt*U
		b = a^2 * (exp(-im*lattice.δt*U) - 1)
		for i in 1:lattice.k-1
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
function sysdynamics_backward!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=DefaultKTruncation)
	# free dynamics
	μ, U = model.μ, model.U
	a = exp(-im*lattice.δt*μ)
	ac = conj(a)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
			pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=band), index(lattice, i+1, conj=false, branch=:-, band=band)
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
				pos1 = index(lattice, i, conj=true, branch=:-, band=band)
				pos2 = index(lattice, i, conj=true, branch=:-, band=band+1)
				pos3 = index(lattice, i+1, conj=false, branch=:-, band=band+1)
				pos4 = index(lattice, i+1, conj=false, branch=:-, band=band)
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=bc)), gmps)
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
sysdynamicsstepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SingleImpurityModel; kwargs...) = si_sysdynamics_stepper!(gmps, lattice; μ=model.μ, U=model.U, kwargs...)

