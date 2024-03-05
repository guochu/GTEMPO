hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::PartialIF; band::Int=1) = hybriddynamics(
				gmps, lattice, corr; band=band, trunc=alg.trunc)
hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::PartialIF; band::Int=1) = hybriddynamics!(
				vacuumstate(lattice), lattice, corr; band=band, trunc=alg.trunc)

hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamics!(GrassmannMPS(scalartype(lattice), length(lattice)), lattice, corr; kwargs...)

"""
	hybriddynamics(gmps::GrassmannMPS, lattice::ImagGrassmannLattice, corr::ImagCorrelationFunction; band, trunc)

imaginary-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	corr = corr1.data
	k = lattice.k-1
	for i in 1:k
		tmp = partialinfluencefunctional(lattice, i+1, [0; view(corr, i, 1:k)], band=band)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps
end

"""
	hybriddynamics(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::MixedCorrelationFunction; band::Int, trunc)

real-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice1Order, corr::MixedCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	kt, Nτ = lattice.kt, lattice.Nτ
	for b1 in (:+, :-)
		for i in 1:kt
			cols_f = [index(corr, i, j, b1=b1, b2=:+) for j in 1:kt]
			cols_b = [index(corr, i, j, b1=b1, b2=:-) for j in 1:kt]
			cols_i = [index(corr, i, j, b1=b1, b2=:τ) for j in 1:Nτ]
			tmp = partialinfluencefunctional(lattice, i, cols_f, cols_b, cols_i, b1=b1, band=band)
			gmps = mult!(gmps, tmp, trunc=trunc)
		end
	end
	b1 = :τ
	for i in 1:Nτ
		cols_f = [index(corr, i, j, b1=b1, b2=:+) for j in 1:kt]
		cols_b = [index(corr, i, j, b1=b1, b2=:-) for j in 1:kt]
		cols_i = [index(corr, i, j, b1=b1, b2=:τ) for j in 1:Nτ]
		tmp = partialinfluencefunctional(lattice, i, cols_f, cols_b, cols_i, b1=b1, band=band)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps		
end


"""
	hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction; band::Int, trunc)

real-time MPS-IF for a single band 
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	k = lattice.k
	for i in 1:k
		tmp1 = partialinfluencefunctional(lattice, i, view(η⁺⁺, i, 1:k), view(η⁺⁻, i, 1:k), b1=:+, band=band)
		tmp3 = partialinfluencefunctional(lattice, i, view(η⁻⁺, i, 1:k), view(η⁻⁻, i, 1:k), b1=:-, band=band)

		gmps = mult!(gmps, tmp1, trunc=trunc)
		gmps = mult!(gmps, tmp3, trunc=trunc)
	end
	return gmps		
end

function hybriddynamicsstepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	@assert length(lattice) == length(gmps) 
	if lattice.k == 2
		pos1, pos2 = index(lattice, 1, conj=true, branch=:+, band=band), index(lattice, 1, conj=false, branch=:+, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=η⁺⁺[1,1])), gmps)
		pos1, pos2 = index(lattice, 1, conj=true, branch=:+, band=band), index(lattice, 1, conj=false, branch=:-, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=η⁺⁻[1,1])), gmps)
		pos1, pos2 = index(lattice, 1, conj=true, branch=:-, band=band), index(lattice, 1, conj=false, branch=:+, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=η⁻⁺[1,1])), gmps)
		pos1, pos2 = index(lattice, 1, conj=true, branch=:-, band=band), index(lattice, 1, conj=false, branch=:-, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=η⁻⁻[1,1])), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
	end
	tmp1 = finalinfluencefunctional(lattice, η⁺⁺, trunc=trunc, b1=:+, b2=:+, band=band)
	tmp2 = finalinfluencefunctional(lattice, η⁺⁻, trunc=trunc, b1=:+, b2=:-, band=band)
	tmp3 = finalinfluencefunctional(lattice, η⁻⁺, trunc=trunc, b1=:-, b2=:+, band=band)
	tmp4 = finalinfluencefunctional(lattice, η⁻⁻, trunc=trunc, b1=:-, b2=:-, band=band)

	gmps = mult!(gmps, tmp1, trunc=trunc)
	gmps = mult!(gmps, tmp2, trunc=trunc)
	gmps = mult!(gmps, tmp3, trunc=trunc)
	gmps = mult!(gmps, tmp4, trunc=trunc)
	return gmps		
end

function hybriddynamicsstepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction; 
		band::Int=1, finalize::Bool=false, trunc::TruncationScheme=DefaultITruncation)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	@assert length(lattice) == length(gmps)
	if lattice.k == 2
		pos1, pos2 = index(lattice, 1, conj=true, branch=:+, band=band), index(lattice, 1, conj=false, branch=:+, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=η⁺⁺[1,1])), gmps)
		pos1, pos2 = index(lattice, 1, conj=true, branch=:+, band=band), index(lattice, 1, conj=false, branch=:-, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=η⁺⁻[1,1])), gmps)
		pos1, pos2 = index(lattice, 1, conj=true, branch=:-, band=band), index(lattice, 1, conj=false, branch=:+, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=η⁻⁺[1,1])), gmps)
		pos1, pos2 = index(lattice, 1, conj=true, branch=:-, band=band), index(lattice, 1, conj=false, branch=:-, band=band)
		apply!(exp(GTerm(pos1, pos2, coeff=η⁻⁻[1,1])), gmps)
		canonicalize!(gmps, alg=Orthogonalize(SVD(), trunc))
	end
	tmp1 = finalinfluencefunctional(lattice, η⁺⁺, finalize=finalize, trunc=trunc, b1=:+, b2=:+, band=band)
	tmp2 = finalinfluencefunctional(lattice, η⁺⁻, finalize=finalize, trunc=trunc, b1=:+, b2=:-, band=band)
	tmp3 = finalinfluencefunctional(lattice, η⁻⁺, finalize=finalize, trunc=trunc, b1=:-, b2=:+, band=band)
	tmp4 = finalinfluencefunctional(lattice, η⁻⁻, finalize=finalize, trunc=trunc, b1=:-, b2=:-, band=band)

	gmps = mult!(gmps, tmp1, trunc=trunc)
	gmps = mult!(gmps, tmp2, trunc=trunc)
	gmps = mult!(gmps, tmp3, trunc=trunc)
	gmps = mult!(gmps, tmp4, trunc=trunc)

	return gmps
end
hybriddynamicsstepper(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = hybriddynamicsstepper!(copy(gmps), lattice, corr; kwargs...)

### for single impurity models

"""
	qim_hybriddynamics(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::NTuple{4, <:AbstractMatrix}; trunc)

Real time hybrid dynamics for single impurity model, all the bands 
share the same bath 
"""
function qim_hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...)
	for band in 1:lattice.bands
		gmps = hybriddynamics!(gmps, lattice, corr; band=band, kwargs...)
	end
	return gmps
end
qim_hybriddynamics(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamics!(copy(gmps), lattice, corr; kwargs...)
qim_hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamics!(vacuumstate(lattice), lattice, corr; kwargs...)

function qim_hybriddynamicsstepper!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...)
	for band in 1:lattice.bands
		gmps = hybriddynamicsstepper!(gmps, lattice, corr; band=band, kwargs...)
	end
	return gmps
end
qim_hybriddynamicsstepper(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction; kwargs...) = qim_hybriddynamicsstepper!(copy(gmps), lattice, corr; kwargs...)

"""
	finalinfluencefunctional(lattice::RealGrassmannLattice2Order, η::AbstractMatrix; kwargs...)

The top right of the influence functional
"""
function finalinfluencefunctional(lattice::RealGrassmannLattice2Order, η::AbstractMatrix; 
				finalize::Bool, trunc::TruncationScheme=DefaultITruncation, kwargs...)
	rows, cols = finalize ? top_right_if_final(lattice, η) : top_right_if(lattice, η)
	out1 = partialinfluencefunctional(lattice, lattice.k, cols; kwargs...)
	out2 = partialinfluencefunctional(lattice, rows, lattice.k; kwargs...)
	return mult(out1, out2, trunc=trunc)
end
function finalinfluencefunctional(lattice::RealGrassmannLattice1Order, η::AbstractMatrix; trunc::TruncationScheme=DefaultITruncation, kwargs...)
	@assert lattice.k <= size(η, 1)
	k = lattice.k
	out1 = partialinfluencefunctional(lattice, k, [η[k, 1:k-1]; 0.5*η[k,k]]; kwargs...)
	out2 = partialinfluencefunctional(lattice, [η[1:k-1, k]; 0.5*η[k,k]], k; kwargs...)
	return mult(out1, out2, trunc=trunc)
end

function top_right_if_final(lattice::RealGrassmannLattice2Order, η::AbstractMatrix)
	@assert lattice.N <= div(size(η,1)-1, 2) 
	k = lattice.k
	cols = zeros(eltype(η), k)
	rows = zeros(eltype(η), k)
	cols[1] = η[2*k-2, 1]
	rows[1] = η[1, 2*k-2]
	for i in 2:k-1
		cols[i] = η[2*k-2, 2*i-2] + η[2*k-2, 2*i-1]
		rows[i] = η[2*i-2, 2*k-2] + η[2*i-1, 2*k-2]
	end
	cols[k] = rows[k] = 0.5 * η[2*k-2, 2*k-2]
	return rows, cols
end

function top_right_if(lattice::RealGrassmannLattice2Order, η::AbstractMatrix)
	@assert lattice.N <= div(size(η,1)-1, 2)
	k = lattice.k
	cols = zeros(eltype(η), k)
	rows = zeros(eltype(η), k)
	cols[1]	= η[2*k-2, 1] + η[2*k-1, 1]
	rows[1] = η[1, 2*k-2] + η[1, 2*k-1]
	for i in 2:k-1
		cols[i] = η[2*k-2, 2*i-2] + η[2*k-2, 2*i-1] + η[2*k-1, 2*i-2] + η[2*k-1, 2*i-1]
		rows[i] = η[2*i-2, 2*k-2] + η[2*i-2, 2*k-1] + η[2*i-1, 2*k-2] + η[2*i-1, 2*k-1]
	end
	tmp = η[2*k-2, 2*k-2] + η[2*k-2, 2*k-1] + η[2*k-1, 2*k-2] + η[2*k-1, 2*k-1]
	cols[k] = rows[k] = 0.5 * tmp
	return rows, cols
end
