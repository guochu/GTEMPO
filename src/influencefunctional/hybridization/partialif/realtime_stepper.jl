function hybriddynamicsstepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	@assert length(lattice) == length(gmps) 
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
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
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
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



"""
	finalinfluencefunctional(lattice::RealGrassmannLattice2Order, η::AbstractMatrix; kwargs...)

The top right of the influence functional
"""
function finalinfluencefunctional(lattice::RealGrassmannLattice2Order, η::AbstractMatrix; 
				finalize::Bool, trunc::TruncationScheme=DefaultITruncation, kwargs...)
	rows, cols = finalize ? top_right_if_final(lattice, η) : top_right_if(lattice, η)
	out1 = partialif_hybrid(lattice, lattice.k, cols; kwargs...)
	out2 = partialif_hybrid(lattice, rows, lattice.k; kwargs...)
	return mult(out1, out2, trunc=trunc)
end
function finalinfluencefunctional(lattice::RealGrassmannLattice1Order, η::AbstractMatrix; trunc::TruncationScheme=DefaultITruncation, kwargs...)
	@assert lattice.k <= size(η, 1)
	k = lattice.k
	out1 = partialif_hybrid(lattice, k, [η[k, 1:k-1]; 0.5*η[k,k]]; kwargs...)
	out2 = partialif_hybrid(lattice, [η[1:k-1, k]; 0.5*η[k,k]], k; kwargs...)
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
