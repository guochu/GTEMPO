"""
	struct AndersonIM

Single-orbital Anderson impurity model with one bath:
H = ϵ_d (n₁ + n₂) + U n₁ n₂ on two bands (one per spin direction).
"""
struct AndersonIM <: ConstImpurityHamiltonian
	U::Float64
	ϵ_d::Float64
end
AndersonIM(U::Real, ϵ_d::Real) = AndersonIM(convert(Float64, U), convert(Float64, ϵ_d))
AndersonIM(; U::Real, ϵ_d::Real) = AndersonIM(U, ϵ_d)

num_bands(::AndersonIM) = 2

"""
	struct ToulouseIM

Toulouse impurity model: the noninteracting (`U = 0`) single-band case of
`AndersonIM`, H = ϵ_d n̂ on a single band.
"""
struct ToulouseIM <: ConstImpurityHamiltonian
	ϵ_d::Float64
end
ToulouseIM(ϵ_d::Real) = ToulouseIM(convert(Float64, ϵ_d))
ToulouseIM(; ϵ_d::Real) = ToulouseIM(ϵ_d)

num_bands(::ToulouseIM) = 1

# function hybriddynamics(gmps::GrassmannMPS, lattice::AbstractImagGrassmannLattice, model::AndersonIM; 
# 					corr::Union{Nothing, <:ImagCorrelationFunction}=nothing, trunc::TruncationScheme=Defaults.alg_trunc_if())
# 	@assert lattice.β == model.bath.β
# 	if isnothing(corr)
# 		corr = correlationfunction(model.bath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end

# function hybriddynamics(gmps::GrassmannMPS, lattice::AbstractRealGrassmannLattice, model::AndersonIM; 
# 					corr::Union{Nothing, <:RealCorrelationFunction}=nothing, 
# 					trunc::TruncationScheme=Defaults.alg_trunc_if())
# 	if isnothing(corr)
# 		corr = correlationfunction(model.bath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end

# function hybriddynamics(gmps::GrassmannMPS, lattice::AbstractMixedGrassmannLattice, model::AndersonIM; 
# 					corr::Union{Nothing, <:MixedCorrelationFunction}=nothing, 
# 					trunc::TruncationScheme=Defaults.alg_trunc_if())
# 	@assert lattice.β == model.bath.β
# 	if isnothing(corr)
# 		corr = correlationfunction(model.bath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end

# function hybriddynamics(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, model::SIDB; 
# 					corr::Union{Nothing, <:ImagCorrelationFunction}=nothing, trunc::TruncationScheme=Defaults.alg_trunc_if())
# 	if isnothing(corr)
# 		corr = correlationfunction(model.leftbath, lattice) + correlationfunction(model.rightbath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end

# function hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::SIDB; 
# 					corr::Union{Nothing, <:RealCorrelationFunction}=nothing, 
# 					trunc::TruncationScheme=Defaults.alg_trunc_if())
# 	if isnothing(corr)
# 		corr = correlationfunction(model.leftbath, lattice) + correlationfunction(model.rightbath, lattice)
# 	end
# 	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
# end


# function sysdynamics!(gmps::GrassmannMPS, lattice::AbstractImagGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=Defaults.alg_trunc_sys()) 
# 	ϵ_d, U = model.ϵ_d, model.U
# 	a, b = siam_coeffs(ϵ_d, U, -lattice.δτ)
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
#             pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
#             apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
#             canonicalize!(gmps, alg=Orthogonalize(orth=SVD(), trunc=trunc))			
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
# 				canonicalize!(gmps, alg=Orthogonalize(orth=SVD(), trunc=trunc))
# 			end
# 		end
# 	end
# 	return gmps	
# end


# function sysdynamics_forward!(gmps::GrassmannMPS, lattice::AbstractRealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=Defaults.alg_trunc_sys())
# 	# free dynamics
# 	ϵ_d, U = model.ϵ_d, model.U
# 	a, b = siam_coeffs(ϵ_d, U, -im*lattice.δt) 
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
#             pos1, pos2 = index(lattice, i+1, conj=true, branch=:+, band=band), index(lattice, i, conj=false, branch=:+, band=band)
#             apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
#             canonicalize!(gmps, alg=Orthogonalize(orth=SVD(), trunc=trunc))			
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
# 				canonicalize!(gmps, alg=Orthogonalize(orth=SVD(), trunc=trunc))			
# 			end
# 		end
# 	end
# 	return gmps
# end
# function sysdynamics_backward!(gmps::GrassmannMPS, lattice::AbstractRealGrassmannLattice, model::SingleImpurityModel; trunc::TruncationScheme=Defaults.alg_trunc_sys())
# 	# free dynamics
# 	ϵ_d, U = model.ϵ_d, model.U
# 	# a = exp(-im*lattice.δt*ϵ_d)
# 	# ac = conj(a)
# 	a, b = siam_coeffs(ϵ_d, U, im*lattice.δt) 
# 	for band in 1:lattice.bands
# 		for i in 1:lattice.k-1
# 			pos1, pos2 = lattice[GrassmannContourIndex(i, :-, band, true)], lattice[GrassmannContourIndex(i+1, :-, band, false)]
# 			apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
# 			canonicalize!(gmps, alg=Orthogonalize(orth=SVD(), trunc=trunc))	
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
# 				canonicalize!(gmps, alg=Orthogonalize(orth=SVD(), trunc=trunc))			
# 			end
# 		end	
# 	end
# 	return gmps
# end

"""
	siam_coeffs(ϵ_d, U, dt) -> (a, b)

Pure-number coefficients of the per-step GTerms of the legacy
Trotterized construction (kept as a numeric utility).
"""
function siam_coeffs(ϵ_d, U, dt)
	# a = exp(dt*(ϵ_d - U/2))
	a = exp(dt*ϵ_d)
	b = a^2 * (exp(dt*U) - 1)
	return a, b
end
