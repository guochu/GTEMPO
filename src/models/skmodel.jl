# The SK impurity model

struct SKIM{B <: AbstractFermionicBath} <: AbstractImpurityModel
	bath::B
	norb::Int
	U::Float64
	J::Float64
	μ::Float64
end
SKIM(bath::AbstractFermionicBath; U::Real, J::Real, norb::Int, μ::Real=-U/2) = SKIM(bath, norb, convert(Float64, U), convert(Float64, J), convert(Float64, μ))


function sysdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, model::SKIM; trunc::TruncationScheme=DefaultKTruncation)
	@assert 2 * model.norb == lattice.bands
	# @assert lattice.β == model.bath.β

	μ, U, J, nimp = model.μ, model.U, model.J, model.norb
	### free dynamics
	a = exp(-lattice.δτ*μ)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
            pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))			
		end
	end

	### interacting dynamics
	a = -lattice.δτ*U
	for x in 1:nimp
		for i in 1:lattice.k - 1
			pos1 = index(lattice, i+1, conj=true, band=2*x)
			pos2 = index(lattice, i+1, conj=true, band=2*x-1)
			pos3 = index(lattice, i, conj=false, band=2*x-1)
			pos4 = index(lattice, i, conj=false, band=2*x)
			apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
			canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))			
		end
	end

	# b = a^2 * (exp(-lattice.δτ*U) - 1)
	# for x in 1:nimp
	# 	for i in 1:lattice.k - 1
	# 		pos1 = index(lattice, i+1, conj=true, band=2*x)
	# 		pos2 = index(lattice, i+1, conj=true, band=2*x-1)
	# 		pos3 = index(lattice, i, conj=false, band=2*x-1)
	# 		pos4 = index(lattice, i, conj=false, band=2*x)
	# 		apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=b)), gmps)
	# 		canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))			
	# 	end
	# end

	a = -lattice.δτ * (U-2*J)
	for x in 1:nimp
		for y in 1:nimp
			if x != y
				for i in 1:lattice.k-1
					pos1 = index(lattice, i+1, conj=true, band=2 * y)
					pos2 = index(lattice, i+1, conj=true, band=2 * x - 1)
					pos3 = index(lattice, i, conj=false, band=2 * x - 1)
					pos4 = index(lattice, i, conj=false, band=2 * y)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
				end
			end
		end
	end

	a = -lattice.δτ * (U-3*J)
	for y in 1:nimp
		for x in (y+1):nimp
			for i in 1:lattice.k-1
				pos1 = index(lattice, i+1, conj=true, band=2*y - 1)
				pos2 = index(lattice, i+1, conj=true, band=2*x - 1)
				pos3 = index(lattice, i, conj=false, band=2*x - 1)
				pos4 = index(lattice, i, conj=false, band=2*y - 1)	
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
			end
			for i in 1:lattice.k-1
				pos1 = index(lattice, i+1, conj=true, band=2*y )
				pos2 = index(lattice, i+1, conj=true, band=2*x )
				pos3 = index(lattice, i, conj=false, band=2*x )
				pos4 = index(lattice, i, conj=false, band=2*y )	
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
			end
		end
	end

	a = lattice.δτ * J
	for x in 1:nimp
		for y in 1:nimp
			if x != y
				for i in 1:lattice.k-1
					pos1 = index(lattice, i+1, conj=true, band=2*x-1)
					pos2 = index(lattice, i+1, conj=true, band=2*x)
					pos3 = index(lattice, i, conj=false, band=2*y-1)
					pos4 = index(lattice, i, conj=false, band=2*y)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
				end
				for i in 1:lattice.k-1
					pos1 = index(lattice, i+1, conj=true, band=2*x-1)
					pos2 = index(lattice, i+1, conj=true, band=2*y)
					pos3 = index(lattice, i, conj=false, band=2*y-1)
					pos4 = index(lattice, i, conj=false, band=2*x)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))					
				end
			end
		end
	end
	return gmps
end


function sysdynamics_forward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, model::SKIM; trunc::TruncationScheme=DefaultKTruncation)
	@assert 2 * model.norb == lattice.bands
	# @assert lattice.β == model.bath.β

	μ, U, J, nimp = model.μ, model.U, model.J, model.norb
	### free dynamics
	a = exp(-im*lattice.δt*μ)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
            pos1, pos2 = index(lattice, i+1, conj=true, band=band, branch=:+), index(lattice, i, conj=false, band=band, branch=:+)
            apply!(exp(GTerm(pos1, pos2, coeff=a)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))			
		end
	end

	### interacting dynamics
	a = -im*lattice.δt*U
	for x in 1:nimp
		for i in 1:lattice.k - 1
			pos1 = index(lattice, i+1, conj=true, band=2*x, branch=:+)
			pos2 = index(lattice, i+1, conj=true, band=2*x-1, branch=:+)
			pos3 = index(lattice, i, conj=false, band=2*x-1, branch=:+)
			pos4 = index(lattice, i, conj=false, band=2*x, branch=:+)
			apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
			canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))			
		end
	end

	a = -im*lattice.δt * (U-2*J)
	for x in 1:nimp
		for y in 1:nimp
			if x != y
				for i in 1:lattice.k-1
					pos1 = index(lattice, i+1, conj=true, band=2 * y, branch=:+)
					pos2 = index(lattice, i+1, conj=true, band=2 * x - 1, branch=:+)
					pos3 = index(lattice, i, conj=false, band=2 * x - 1, branch=:+)
					pos4 = index(lattice, i, conj=false, band=2 * y, branch=:+)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
				end
			end
		end
	end

	a = -im*lattice.δt * (U-3*J)
	for y in 1:nimp
		for x in (y+1):nimp
			for i in 1:lattice.k-1
				pos1 = index(lattice, i+1, conj=true, band=2*y - 1, branch=:+)
				pos2 = index(lattice, i+1, conj=true, band=2*x - 1, branch=:+)
				pos3 = index(lattice, i, conj=false, band=2*x - 1, branch=:+)
				pos4 = index(lattice, i, conj=false, band=2*y - 1, branch=:+)	
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
			end
			for i in 1:lattice.k-1
				pos1 = index(lattice, i+1, conj=true, band=2*y, branch=:+)
				pos2 = index(lattice, i+1, conj=true, band=2*x, branch=:+)
				pos3 = index(lattice, i, conj=false, band=2*x, branch=:+)
				pos4 = index(lattice, i, conj=false, band=2*y, branch=:+)	
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
			end
		end
	end

	a = im*lattice.δt * J
	for x in 1:nimp
		for y in 1:nimp
			if x != y
				for i in 1:lattice.k-1
					pos1 = index(lattice, i+1, conj=true, band=2*x-1, branch=:+)
					pos2 = index(lattice, i+1, conj=true, band=2*x, branch=:+)
					pos3 = index(lattice, i, conj=false, band=2*y-1, branch=:+)
					pos4 = index(lattice, i, conj=false, band=2*y, branch=:+)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
				end
				for i in 1:lattice.k-1
					pos1 = index(lattice, i+1, conj=true, band=2*x-1, branch=:+)
					pos2 = index(lattice, i+1, conj=true, band=2*y, branch=:+)
					pos3 = index(lattice, i, conj=false, band=2*y-1, branch=:+)
					pos4 = index(lattice, i, conj=false, band=2*x, branch=:+)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=a)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))					
				end
			end
		end
	end

	return gmps
end

function sysdynamics_backward!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, model::SKIM; trunc::TruncationScheme=DefaultKTruncation)
	@assert 2 * model.norb == lattice.bands
	# @assert lattice.β == model.bath.β

	μ, U, J, nimp = model.μ, model.U, model.J, model.norb
	### free dynamics
	a = exp(-im*lattice.δt*μ)
	ac = conj(a)
	for band in 1:lattice.bands
		for i in 1:lattice.k-1
            pos1, pos2 = index(lattice, i, conj=true, band=band, branch=:-), index(lattice, i+1, conj=false, band=band, branch=:-)
            apply!(exp(GTerm(pos1, pos2, coeff=ac)), gmps)
            canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))			
		end
	end	

	### interacting dynamics
	a = -im*lattice.δt*U
	ac = conj(a)
	for x in 1:nimp
		for i in 1:lattice.k - 1
			pos1 = index(lattice, i, conj=true, band=2*x, branch=:-)
			pos2 = index(lattice, i, conj=true, band=2*x-1, branch=:-)
			pos3 = index(lattice, i+1, conj=false, band=2*x-1, branch=:-)
			pos4 = index(lattice, i+1, conj=false, band=2*x, branch=:-)
			apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=ac)), gmps)
			canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))			
		end
	end	


	a = -im*lattice.δt * (U-2*J)
	ac = conj(a)
	for x in 1:nimp
		for y in 1:nimp
			if x != y
				for i in 1:lattice.k-1
					pos1 = index(lattice, i, conj=true, band=2 * y, branch=:-)
					pos2 = index(lattice, i, conj=true, band=2 * x - 1, branch=:-)
					pos3 = index(lattice, i+1, conj=false, band=2 * x - 1, branch=:-)
					pos4 = index(lattice, i+1, conj=false, band=2 * y, branch=:-)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=ac)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
				end
			end
		end
	end

	a = -im*lattice.δt * (U-3*J)
	ac = conj(a)
	for y in 1:nimp
		for x in (y+1):nimp
			for i in 1:lattice.k-1
				pos1 = index(lattice, i, conj=true, band=2*y - 1, branch=:-)
				pos2 = index(lattice, i, conj=true, band=2*x - 1, branch=:-)
				pos3 = index(lattice, i+1, conj=false, band=2*x - 1, branch=:-)
				pos4 = index(lattice, i+1, conj=false, band=2*y - 1, branch=:-)	
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=ac)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
			end
			for i in 1:lattice.k-1
				pos1 = index(lattice, i, conj=true, band=2*y, branch=:-)
				pos2 = index(lattice, i, conj=true, band=2*x, branch=:-)
				pos3 = index(lattice, i+1, conj=false, band=2*x, branch=:-)
				pos4 = index(lattice, i+1, conj=false, band=2*y, branch=:-)	
				apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=ac)), gmps)
				canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
			end
		end
	end

	a = im*lattice.δt * J
	ac = conj(a)
	for x in 1:nimp
		for y in 1:nimp
			if x != y
				for i in 1:lattice.k-1
					pos1 = index(lattice, i, conj=true, band=2*x-1, branch=:-)
					pos2 = index(lattice, i, conj=true, band=2*x, branch=:-)
					pos3 = index(lattice, i+1, conj=false, band=2*y-1, branch=:-)
					pos4 = index(lattice, i+1, conj=false, band=2*y, branch=:-)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=ac)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))
				end
				for i in 1:lattice.k-1
					pos1 = index(lattice, i, conj=true, band=2*x-1, branch=:-)
					pos2 = index(lattice, i, conj=true, band=2*y, branch=:-)
					pos3 = index(lattice, i+1, conj=false, band=2*y-1, branch=:-)
					pos4 = index(lattice, i+1, conj=false, band=2*x, branch=:-)	
					apply!(exp(GTerm(pos1, pos2, pos3, pos4, coeff=ac)), gmps)
					canonicalize!(gmps, alg=Orthogonalize(trunc = trunc))					
				end
			end
		end
	end
	return gmps
end




function hybriddynamics(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, model::SKIM; 
					corr::Union{Nothing, <:ImagCorrelationFunction}=nothing, trunc::TruncationScheme=DefaultITruncation)
	@assert 2 * model.norb == lattice.bands
	@assert lattice.β == model.bath.β
	if isnothing(corr)
		corr = correlationfunction(model.bath, lattice)
	end
	return qim_hybriddynamics(gmps, lattice, corr, trunc=trunc)
end

