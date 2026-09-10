

"""
    hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::ExactTTIIF; band::Int=1)

Construct the influence functional with the `ExactTTIIF` algorithm on
`lattice` (starting from the vacuum state). Equivalently the in-place
[`hybriddynamics!`](@ref), see its docstring for the workflow of merging the
influence functional into an existing `GrassmannMPS` such as the impurity
dynamics from `sysdynamics`.
"""
hybriddynamics(lattice::AbstractGrassmannLattice, 
corr::AbstractCorrelationFunction, alg::ExactTTIIF; kwargs...) = hybriddynamics!(vacuumstate(lattice), lattice, corr, alg; kwargs...)

"""
    hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::ExactTTIIF; band::Int=1)

In-place version of the `ExactTTIIF` algorithm: the influence functional is
multiplied onto `gmps` directly (term by term), so that `gmps` carries both
the impurity dynamics and the bath influence. A typical workflow is

    gmps = sysdynamics(lat, model, trunc=trunc)     # impurity dynamics
    hybriddynamics!(gmps, lat, corr, ExactTTIIF())  # merge the IF in place

For `lattice.bands > 1` the influence functional terms are built on a
single-band lattice and expanded to the given `band` via `fillband`.
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::ExactTTIIF; band::Int=1)
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	mpss = _influencefunctional(lattice, corr, alg; band=band)
	for i in 1:length(mpss)
		t = @elapsed gmps = mult!(gmps, mpss[i], alg.algmult)
		(alg.verbosity >= 2) && println("$i of $(length(mpss)) takes $t seconds, result mps of bond dimension: ", bond_dimension(gmps))
	end
	return gmps
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





