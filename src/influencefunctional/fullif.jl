# Full influence functional

function hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::TranslationInvariantIF; 
						band::Int=1, trunc::TruncationDimCutoff=DefaultITruncation)
	mps = hybriddynamics(lattice, corr, alg, band=band, trunc=trunc)
	return mult!(gmps, mps, trunc=trunc)
end

function hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::TranslationInvariantIF; 
						band::Int=1, trunc::TruncationDimCutoff=DefaultITruncation)
	trunc0 = truncdimcutoff(D=trunc.D, ϵ=trunc.ϵ/2^(alg.k), add_back=0)
	mps = differentialinfluencefunctional(lattice, corr, 1/2^(alg.k), alg.algevo, band=band, algexpan=alg.algexpan, trunc=trunc0)
	return _exp_mult(mps, alg.k, trunc)
end

function _exp_mult(mps, k::Int, trunc::TruncationDimCutoff)
	D, ϵ0 = trunc.D, trunc.ϵ
	scale = 2
	dtt = 1/scale^k	
	for i in 1:k
		dtt *= scale
		mps = mult(mps, mps, trunc=truncdimcutoff(D=D, ϵ=ϵ0*dtt, add_back=0))
		# println("mps bond dimension is ", bond_dimension(mps), " at ", i, "-th iteration")
	end	
	return mps
end

