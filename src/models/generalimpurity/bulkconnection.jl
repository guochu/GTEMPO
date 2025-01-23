"""
	bulkconnection(x::GrassmannMPS, lattice::AbstractGrassmannLattice; kwargs...)

Similar to boundarycondition, but applies the bulk connection terms instead of 
the boundary connection term
"""
bulkconnection(lattice::AbstractGrassmannLattice; kwargs...) = bulkconnection!(vacuumstate(lattice), lattice; kwargs...)
# imaginary-time
function bulkconnection!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	return bulkconnection_util!(gmps, lattice, lattice.Nτ, band, :τ, trunc)
end
function bulkconnection!(gmps::GrassmannMPS, lattice::RealGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	gmps = bulkconnection_util!(gmps, lattice, lattice.Nt, band, :+, trunc)
	return bulkconnection_util!(gmps, lattice, lattice.Nt, band, :-, trunc)
end
function bulkconnection!(gmps::GrassmannMPS, lattice::MixedGrassmannLattice; band::Int=1, trunc::TruncationScheme=DefaultIntegrationTruncation)
	gmps = bulkconnection_util!(gmps, lattice, lattice.Nt, band, :+, trunc)
	gmps = bulkconnection_util!(gmps, lattice, lattice.Nt, band, :-, trunc)
	return bulkconnection_util!(gmps, lattice, lattice.Nτ, band, :τ, trunc)
end


# the overlap between ⟨āⱼ₊₁ āⱼ
function bulkconnection_util!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, N::Int, band::Int, branch::Symbol, trunc)
	# applying the overlap of GVs
	alg = Orthogonalize(trunc = trunc)
	for j in 1:N
		a, b = (branch == :-) ? (j, j+1) : (j+1, j)
		pos1, pos2 = index(lattice, a, conj=true, band=band, branch=branch), index(lattice, b, conj=false, band=band, branch=branch)
		m = exp(GTerm(pos1, pos2, coeff=1))
		apply!(m, gmps)
		canonicalize!(gmps, alg=alg)
	end
	return gmps
end