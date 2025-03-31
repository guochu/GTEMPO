
# naive version
function hybriddynamics_naive!(gmps::FockMPS, lattice::MixedFockLattice1Order, corr::AbstractMixedCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	(lattice.bands in (1, 2)) || throw(ArgumentError("number of bands should be either 1 or 2"))
	if lattice.bands == 1
		return hybriddynamics_1band_naive!(gmps, lattice, corr, trunc=trunc)
	else
		return hybriddynamics_2band_naive!(gmps, lattice, corr, trunc=trunc)
	end
end

function hybriddynamics_1band_naive!(gmps::FockMPS, lattice::MixedFockLattice1Order, corr::AbstractMixedCorrelationFunction; kwargs...)
	@assert lattice.bands == 1
	return _hybriddynamics_1band_naive!(gmps, lattice, corr, 1; kwargs...)
end 


function _hybriddynamics_1band_naive!(gmps::FockMPS, lattice::MixedFockLattice1Order, corr::AbstractMixedCorrelationFunction, band::Int; 
												trunc::TruncationScheme=DefaultITruncation)	
	# (LayoutStyle(lattice) isa BandLocalLayout) || throw(ArgumentError("currently only TimelocalLayout support for this function"))
	alg = Orthogonalize(TK.SVD(), trunc)
	for b1 in branches(lattice), b2 in branches(lattice)
		k1 = ifelse(b1==:τ, lattice.Nτ, lattice.Nt)
		k2 = ifelse(b2==:τ, lattice.Nτ, lattice.Nt)
		for i in 1:k1, j in 1:k2
			pos1, pos2 = index(lattice, i, band=band, branch=b1), index(lattice, j, band=band, branch=b2)
			c = exp(index(corr, i, j, b1=b1, b2=b2)) - 1
			if pos1 == pos2
				t = exp(NTerm(pos1, coeff=c)) 
			else
				t = exp(NTerm(pos1, pos2, coeff=c)) 
			end
			apply!(t, gmps)
			canonicalize!(gmps, alg=alg)			
		end
	end
	return gmps
end

function hybriddynamics_2band_naive!(gmps::FockMPS, lattice::MixedFockLattice1Order, corr::AbstractMixedCorrelationFunction; 
												trunc::TruncationScheme=DefaultMPOTruncation)
	@assert lattice.bands == 2
	for band in 1:lattice.bands
		_hybriddynamics_1band_naive!(gmps, lattice, corr, band, trunc=trunc)
	end
	alg = Orthogonalize(TK.SVD(), trunc)

	for b1 in branches(lattice), b2 in branches(lattice)
		k1 = ifelse(b1==:τ, lattice.Nτ, lattice.Nt)
		k2 = ifelse(b2==:τ, lattice.Nτ, lattice.Nt)
		for i in 1:k1, j in 1:k2
			pos1, pos2 = index(lattice, i, band=1, branch=b1), index(lattice, j, band=2, branch=b2)
			c = exp(2*index(corr, i, j, b1=b1, b2=b2)) - 1
			if pos1 == pos2
				t = exp(NTerm(pos1, coeff=c)) 
			else
				t = exp(NTerm(pos1, pos2, coeff=c)) 
			end
			apply!(t, gmps)
			canonicalize!(gmps, alg=alg)						
		end
	end	

	return gmps
end
