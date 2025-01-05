function sysdynamics_imaginary!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
	(lattice.bands == h.bands) || throw(ArgumentError("lattice bands and ImpurityHamiltonian bands mismatch"))
	return sysdynamics_util!(gmps, lattice, h, lattice.N, -lattice.δτ, :τ, trunc)
end
function sysdynamics_forward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
	(lattice.bands == h.bands) || throw(ArgumentError("lattice bands and ImpurityHamiltonian bands mismatch"))
	return sysdynamics_util!(gmps, lattice, h, lattice.Nt, -im*lattice.δt, :+, trunc)
end
function sysdynamics_backward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
	(lattice.bands == h.bands) || throw(ArgumentError("lattice bands and ImpurityHamiltonian bands mismatch"))
	return sysdynamics_util!(gmps, lattice, h, lattice.Nt, im*lattice.δt, :-, trunc)
end

baresysdynamics(lattice::AbstractGrassmannLattice, h::ImpurityHamiltonian; kwargs...) = baresysdynamics!(vacuumstate(lattice), lattice, h; kwargs...)

"""
	sysdynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::ImpurityHamiltonian; trunc)
"""
function baresysdynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, h::ImpurityHamiltonian; trunc::TruncationScheme=DefaultKTruncation)
	(lattice.bands == h.bands) || throw(ArgumentError("lattice bands and ImpurityHamiltonian bands mismatch"))
	return baresysdynamics_util!(gmps, lattice, h, lattice.N, -lattice.δτ, :τ, trunc)
end

function baresysdynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, h::ImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	(lattice.bands == h.bands) || throw(ArgumentError("lattice bands and ImpurityHamiltonian bands mismatch"))
	if isnothing(branch)
		baresysdynamics_util!(gmps, lattice, h, lattice.N, -im*lattice.δt, :+, trunc)
		return baresysdynamics_util!(gmps, lattice, h, lattice.N, im*lattice.δt, :-, trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		if branch == :+
			return baresysdynamics_util!(gmps, lattice, h, lattice.N, -im*lattice.δt, :+, trunc)
		else
			return baresysdynamics_util!(gmps, lattice, h, lattice.N, im*lattice.δt, :-, trunc)
		end
	end
end

function baresysdynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, h::ImpurityHamiltonian; 
						branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	(lattice.bands == h.bands) || throw(ArgumentError("lattice bands and ImpurityHamiltonian bands mismatch"))
	if isnothing(branch)
		baresysdynamics_util!(gmps, lattice, h, lattice.Nt, -im*lattice.δt, :+, trunc)
		baresysdynamics_util!(gmps, lattice, h, lattice.Nt, im*lattice.δt, :-, trunc)
		return baresysdynamics_util!(gmps, lattice, h, lattice.Nτ, -lattice.δτ, :τ, trunc)
	else
		(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
		if branch == :+
			return baresysdynamics_util!(gmps, lattice, h, lattice.Nt, -im*lattice.δt, :+, trunc)
		elseif branch == :-
			return baresysdynamics_util!(gmps, lattice, h, lattice.Nt, im*lattice.δt, :-, trunc)
		else
			return baresysdynamics_util!(gmps, lattice, h, lattice.Nτ, -lattice.δτ, :τ, trunc)
		end
	end
end


function sysdynamics_util!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::ImpurityHamiltonian, N::Int, dt::Number, branch::Symbol, trunc)
	for j in 1:N
		baresysdynamics_onestep_util!(gmps, lattice, h, j, dt, branch, trunc)
	end
	for band in 1:lattice.bands
		bulkconnection_util!(gmps, lattice, N, band, branch, trunc)
	end
	return gmps
end

function baresysdynamics_util!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::ImpurityHamiltonian, N::Int, dt::Number, branch::Symbol, trunc)
	for j in 1:N
		baresysdynamics_onestep_util!(gmps, lattice, h, j, dt, branch, trunc)
	end
	return gmps
end


function baresysdynamics_onestep_util!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, h::ImpurityHamiltonian, j::Int, dt::Number, branch::Symbol, trunc)
	alg = Orthogonalize(trunc = trunc)
	for hj in h.data
		t = get_Gterm(lattice, hj, j, dt, branch)
		apply!(t, gmps)
		canonicalize!(gmps, alg=alg)		
	end
	return gmps
end

# aᵢ+ aⱼ
function get_Gterm(lattice::AbstractGrassmannLattice, h::TunnelingTerm, j::Int, dt::Number, branch::Symbol)
	b1, b2 = positions(h)
	c = h.coeff * dt
	a, b = (branch == :-) ? (j, j+1) : (j+1, j)
	pos1, pos2 = index(lattice, a, conj=true, band=b1, branch=branch), index(lattice, b, conj=false, band=b2, branch=branch)
	return exp(GTerm(pos1, pos2, coeff=c))
end

function get_Gterm(lattice::AbstractGrassmannLattice, h::InteractionTerm, j::Int, dt::Number, branch::Symbol)
	b1, b2, b3, b4 = positions(h)
	c = h.coeff * dt
	a, b = (branch == :-) ? (j, j+1) : (j+1, j)
	pos1, pos2 = index(lattice, a, conj=true, band=b1, branch=branch), index(lattice, a, conj=true, band=b2, branch=branch)
	pos3, pos4 = index(lattice, b, conj=false, band=b3, branch=branch), index(lattice, b, conj=false, band=b4, branch=branch)
	return exp(GTerm(pos1, pos2, pos3, pos4, coeff=c))
end
