


function hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::ExactTranslationInvariantIF; band::Int=1)
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	if lattice.bands == 1
		return differentialinfluencefunctional(lattice, corr, alg)
	else
		lattice1 = similar(lattice, bands=1)
		mps = differentialinfluencefunctional(lattice1, corr, alg)
		return fillband(lattice, mps; band=band)
	end
end



# get exact WII
function exp_QZ_ZQ(t::SchurMPOTensor)
	λ = t[2,2].data[end]
	a = t[2,3].data[2]
	@assert a == λ
	b = t[1,2].data[5]
	c = t[1,3].data[2]

	M = Matrix{Any}(undef, 2, 2)
	M[1,1] = deepcopy(t[1,3])
	M[1,2] = deepcopy(t[1,2])
	M[2,1] = deepcopy(t[2,3])
	M[2,2] = deepcopy(t[2,2])
	M[1,1].data .+= [1, 0, 0, 1, 1, 0, 0, 1]
	M[2,2].data[6] += (-a*b + c*λ)

	return SparseMPOTensor(M), λ, b, c
end

function exp_IQ_QI(t::SchurMPOTensor)
	λ = t[2,2].data[end]
	a = t[2,3].data[1]
	@assert a == λ
	b = t[1,2].data[2]
	c = t[1,3].data[2]

	M = Matrix{Any}(undef, 2, 2)
	M[1,1] = deepcopy(t[1,3])
	M[1,2] = deepcopy(t[1,2])
	M[2,1] = deepcopy(t[2,3])
	M[2,2] = deepcopy(t[2,2])
	M[1,1].data .+= [1, 0, 0, 1, 1, 0, 0, 1]
	M[2,2].data[6] += (a*b + c*λ)
	return SparseMPOTensor(M), λ, b, c
end

function exact_WII(t::SchurMPOTensor)
	t12 = t[1,2].data
	if t12[2] != 0
		return exp_IQ_QI(t)
	elseif t12[4] != 0
		return exp_QZ_ZQ(t)
	else
		error("Invalid input for exact_WII")
	end
end


include("imaginarytime.jl")
include("realtime.jl")





