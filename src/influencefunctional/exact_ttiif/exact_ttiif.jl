

"""
    hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::ExactTTIIF; band::Int=1)

Construct the influence functional with the `ExactTTIIF` algorithm on
`lattice` (starting from the vacuum state). Equivalently the in-place
`hybriddynamics!`, see its docstring for the workflow of merging the
influence functional into an existing `GrassmannMPS` such as the impurity
dynamics from `sysdynamics`.
"""
function hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::ExactTTIIF; kwargs...)
	T = promote_type(scalartype(lattice), scalartype(corr))
	return hybriddynamics!(vacuumstate(T, lattice), lattice, corr, alg; kwargs...)
end

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
	# the vertex-shifted Prony channels can be complex even for real
	# correlation data: promote the target MPS instead of failing in mult!
	scalartype(gmps) == scalartype(mpss[1]) || (gmps = complex(gmps))
	for i in 1:length(mpss)
		t = @elapsed gmps = mult!(gmps, mpss[i], alg.alg_mult)[1]
		(alg.verbosity >= 2) && println("$i of $(length(mpss)) takes $t seconds, result mps of bond dimension: ", bonddim(gmps))
	end
	# vertex-shifted convention: multiply the turning-point window factors
	for c in (corr isa RealCorrelationFunction ? _window_corrections(lattice, corr, band=band) : ())
		gmps = mult!(gmps, c, alg.alg_mult)[1]
	end
	# vertex-shifted convention: exact single-pair gates on whatever the bulk
	# + window factors leave on the local (equal-node / nearest-neighbour
	# contact) and windowed turning-point positions. The exp of one Grassmann
	# pair is nilpotent (1 + pair), so each gate is exact; the coefficients
	# are measured on a scratch product of the IF factors ALONE (never on the
	# incoming gmps — an impurity dynamics merged in-place carries its own
	# bath-leg amplitudes that are not part of the strict IF kernel). This is
	# the validated repair recipe (the gtempo-shift staged construction) with
	# the measurement in place of the fitted-pole compensation; the extra
	# scratch product doubles the IF-side mult count per call.
	if corr isa RealCorrelationFunction
		scratch = vacuumstate(promote_type(scalartype(mpss[1]), scalartype(corr)), lattice)
		for m in mpss
			scratch = mult!(scratch, m, alg.alg_mult)[1]
		end
		for c in _window_corrections(lattice, corr, band=band)
			scratch = mult!(scratch, c, alg.alg_mult)[1]
		end
		gates = _realtime_residual_gates(scratch, lattice, corr)
		(alg.verbosity >= 1) && println(length(gates), " residual gates to be applied...")
		for g in gates
			gmps = mult!(gmps, g, alg.alg_mult)[1]
		end
	end
	return gmps
end

# --------------------------------------------------------------------------- #
# residual local/boundary gates of the vertex-shifted realtime ExactTTIIF:
# the shifted bulk reproduces the strict node-vertex kernel on the smooth
# Toeplitz positions exactly (Prony tails + native nearest-neighbour contact
# row), while the windowed turning-point legs keep subtraction residuals of
# the `_window_correction` factors. The strict target kernel is read through
# the correlation-function accessor (the same convention the PartialIF
# assembly consumes, cross-checked to float precision), the scratch IF
# product's kernel is read back by contracting the physical blocks, and every
# residual entry is closed by an exact `exp(GTerm)` gate.
# --------------------------------------------------------------------------- #
function _strict_realtime_kernel(lattice::AbstractRealGrassmannLattice, corr::RealCorrelationFunction)
	pairs = Tuple{GrassmannContourIndex,GrassmannContourIndex,scalartype(corr)}[]
	for band in 1:lattice.bands, br1 in (:+, :-), j1 in 1:lattice.k, br2 in (:+, :-), j2 in 1:lattice.k
		r = GrassmannContourIndex(j1, br1, band, true)    # ā leg
		c = GrassmannContourIndex(j2, br2, band, false)   # a leg
		v = corr[conj(r), conj(c)] * (br1 !== br2 ? -1 : 1)
		push!(pairs, (r, c, v))
	end
	return pairs
end

# dense (2×2 per occupation) block of one site tensor, contracted over the
# Z2 fusion trees
function _dense_site_matrix(A, σ::Int)
	Vl, Vr = space_l(A), space_r(A)
	M = zeros(scalartype(A), Z2Tensors.dim(Vl), Z2Tensors.dim(Vr))
	for (fl, fr) in Z2Tensors.fusiontrees(A)
		sl, sp = fl.uncoupled
		sp.n == σ || continue
		fl.coupled.n == fr.coupled.n || continue
		sr = fr.uncoupled[1]
		o1 = sum(d.n < sl.n ? Z2Tensors.dim(Vl, d) : 0 for d in Z2Tensors.sectors(Vl))
		o2 = sum(d.n < sr.n ? Z2Tensors.dim(Vr, d) : 0 for d in Z2Tensors.sectors(Vr))
		dl, dr = Z2Tensors.dim(Vl, sl), Z2Tensors.dim(Vr, sr)
		M[o1+1:o1+dl, o2+1:o2+dr] .= reshape(collect(A[fl, fr]), dl, dr)
	end
	return M
end

# bath-bath bilinear kernel of an assembled GrassmannMPS: read the two-leg
# occupation amplitudes off the contracted site blocks (normalised by the
# vacuum amplitude; the state scaling cancels in the ratio)
function _assembled_kernel(psi::GrassmannMPS, lattice::AbstractRealGrassmannLattice, pairs)
	ms = [[_dense_site_matrix(psi[p], σ) for σ in 0:1] for p in 1:length(psi)]
	coeff(occ) = begin
		v = ones(ComplexF64, 1, 1)
		for p in eachindex(ms)
			v = v * ms[p][p in occ ? 2 : 1]
		end
		return v[1, 1]
	end
	c0 = coeff(Int[])
	D = Vector{ComplexF64}(undef, length(pairs))
	for (n, (r, c, _)) in enumerate(pairs)
		pr, pc = lattice[r], lattice[c]
		D[n] = (pr < pc ? 1 : -1) * coeff([pr, pc]) / c0
	end
	return D
end

function _realtime_residual_gates(psi::GrassmannMPS, lattice::AbstractRealGrassmannLattice, corr::RealCorrelationFunction)
	pairs = _strict_realtime_kernel(lattice, corr)
	D = _assembled_kernel(psi, lattice, pairs)
	scale = maximum(abs(v[3]) for v in pairs)
	tol = 1.0e-12 * max(scale, 1.0)
	gates = GrassmannMPS[]
	for (n, (r, c, v)) in enumerate(pairs)
		Δ = v - D[n]
		abs(Δ) <= tol && continue
		pr, pc = lattice[r], lattice[c]
		g = vacuumstate(promote_type(scalartype(psi), typeof(Δ)), lattice)
		apply!(exp(GTerm(pr, pc; coeff=Δ)), g)
		push!(gates, g)
	end
	return gates
end


# get exact WII
function exp_QZ_ZQ(t::SchurMPOTensor)
	λ = t[2,2].data[end]
	a = t[2,3].data[2]
	# a == λ identifies the pole on the closing operator; λ = 0 is the native
	# nearest-neighbour contact channel (NearestNeighbourTerm): the closing
	# operator carries no pole and the same closed form holds (verified
	# against the WII exponential to float precision)
	@assert (a == λ || iszero(λ))
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
	# see exp_QZ_ZQ: λ = 0 accepts the pole-less nearest-neighbour channel
	@assert (a == λ || iszero(λ))
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





