"""
	hybriddynamicsstepper!(gmps, lattice, corr; band=1, trunc=Defaults.alg_trunc_if())

In-place single time step of the influence functional dynamics: applies the
IF terms connecting the current time step to all previous steps onto `gmps`
and moves it forward by one step. The non-mutating version
`hybriddynamicsstepper` operates on a copy.
"""
function hybriddynamicsstepper!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, corr::RealCorrelationFunction; band::Int=1, trunc::TruncationScheme=Defaults.alg_trunc_if())
	η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
	@assert size(η⁺⁺) == size(η⁺⁻) == size(η⁻⁺) == size(η⁻⁻)
	@assert length(lattice) == length(gmps) 
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	if lattice.k == 2
		pos1, pos2 = lattice[GrassmannContourIndex(1, :+, band, true)], lattice[GrassmannContourIndex(1, :+, band, false)]
		apply!(exp(GTerm(pos1, pos2, coeff=η⁺⁺[1,1])), gmps)
		pos1, pos2 = lattice[GrassmannContourIndex(1, :+, band, true)], lattice[GrassmannContourIndex(1, :-, band, false)]
		apply!(exp(GTerm(pos1, pos2, coeff=η⁺⁻[1,1])), gmps)
		pos1, pos2 = lattice[GrassmannContourIndex(1, :-, band, true)], lattice[GrassmannContourIndex(1, :+, band, false)]
		apply!(exp(GTerm(pos1, pos2, coeff=η⁻⁺[1,1])), gmps)
		pos1, pos2 = lattice[GrassmannContourIndex(1, :-, band, true)], lattice[GrassmannContourIndex(1, :-, band, false)]
		apply!(exp(GTerm(pos1, pos2, coeff=η⁻⁻[1,1])), gmps)
		canonicalize!(gmps, alg=Orthogonalize(orth=SVD(), trunc=trunc))
	end
	tmp1 = finalinfluencefunctional(lattice, η⁺⁺, trunc=trunc, b1=:+, b2=:+, band=band)
	tmp2 = finalinfluencefunctional(lattice, η⁺⁻, trunc=trunc, b1=:+, b2=:-, band=band)
	tmp3 = finalinfluencefunctional(lattice, η⁻⁺, trunc=trunc, b1=:-, b2=:+, band=band)
	tmp4 = finalinfluencefunctional(lattice, η⁻⁻, trunc=trunc, b1=:-, b2=:-, band=band)

	gmps = mult!(gmps, tmp1, trunc=trunc)[1]
	gmps = mult!(gmps, tmp2, trunc=trunc)[1]
	gmps = mult!(gmps, tmp3, trunc=trunc)[1]
	gmps = mult!(gmps, tmp4, trunc=trunc)[1]
	return gmps
end

function finalinfluencefunctional(lattice::RealGrassmannLattice, η::AbstractMatrix; trunc::TruncationScheme=Defaults.alg_trunc_if(), kwargs...)
	@assert lattice.k <= size(η, 1)
	k = lattice.k
	out1 = partialif_hybrid(lattice, k, [η[k, 1:k-1]; 0.5*η[k,k]]; kwargs...)
	out2 = partialif_hybrid(lattice, [η[1:k-1, k]; 0.5*η[k,k]], k; kwargs...)
	return mult(out1, out2, trunc=trunc)[1]
end


