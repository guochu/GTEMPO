# Full influence functional

function hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::TranslationInvariantIF; band::Int=1)
	mps = hybriddynamics(lattice, corr, alg, band=band)
	return mult!(gmps, mps, alg.algmult)
end


function hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::TranslationInvariantIF; band::Int=1)
	algmult0 = alg.algmult
	trunc0 = algmult0.trunc
	D, ϵ0 = trunc0.D, trunc0.ϵ
	algmult = changetrunc(algmult0, trunc=truncdimcutoff(D=D, ϵ=ϵ0/2^(alg.k)))
	# mps = differentialinfluencefunctional(lattice, corr, 1/2^(alg.k), alg.algevo, algmult, band=band, algexpan=alg.algexpan)
	if alg.verbosity > 1
		t = @elapsed mps = differentialinfluencefunctional(lattice, corr, 1/2^(alg.k), alg.algevo, algmult, band=band, algexpan=alg.algexpan)
		println("building the initial MPS-IF takes $t seconds, bond dimension is ", bond_dimension(mps))
	else
		mps = differentialinfluencefunctional(lattice, corr, 1/2^(alg.k), alg.algevo, algmult, band=band, algexpan=alg.algexpan)
	end
	

	scale = 2
	dtt = 1/scale^(alg.k)
	for i in 1:alg.k
		dtt *= scale
		# mps = mult(mps, mps, trunc=truncdimcutoff(D=D, ϵ=ϵ0*dtt, add_back=0))
		trunc = truncdimcutoff(D=D, ϵ=ϵ0 * dtt)
		algmult = changetrunc(algmult0, trunc=trunc)
		# mps = mult(mps, mps, algmult)
		if alg.verbosity > 1
			t = @elapsed mps = mult(mps, mps, algmult)
			println("the $i-th iteration takes $t seconds, bond dimension is ", bond_dimension(mps))
		else
			mps = mult(mps, mps, algmult)
		end
		# println("mps bond dimension is ", bond_dimension(mps), " at ", i, "-th iteration")
	end	
	return mps
end

# function _exp_mult(mps, k::Int, trunc::TruncationDimCutoff)
# 	D, ϵ0 = trunc.D, trunc.ϵ
# 	scale = 2
# 	dtt = 1/scale^k	
# 	for i in 1:k
# 		dtt *= scale
# 		# mps = mult(mps, mps, trunc=truncdimcutoff(D=D, ϵ=ϵ0*dtt, add_back=0))
# 		mps = iterativemult(mps, mps, trunc=truncdimcutoff(D=D, ϵ=ϵ0*dtt, add_back=0))
# 		# println("mps bond dimension is ", bond_dimension(mps), " at ", i, "-th iteration")
# 	end	
# 	return mps
# end

