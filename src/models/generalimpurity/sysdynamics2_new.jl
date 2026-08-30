# sysdynamics2_new
# ---------------
# Same propagator as sysdynamics2, but the Fock-space propagator
# exp(coeff*H) is converted *directly* into the site tensors of a
# GrassmannMPS in the coherent-state (Grassmann) representation,
# instead of being decomposed into GTerms which are then applied
# one by one to the vacuum state.
#
# The coherent-state coefficient tensor
#     T[x] = Σ_{m,n} (-1)^{inv(m,n)} ρ[m,n]
# (x = occupations of the bra/ket Grassmann variables in site order,
# ρ = fockstate, inv = inversions needed to sort the monomial
#  c̄(bra, ascending bands) c(ket, descending bands) into site order)
# is assembled directly in the Z2Irrep(0) block of a graded tensor and
# then factorized into site tensors by a single Z2-graded SVD sweep.
# The propagator of every time step has the same structure, so the site
# tensors are cached per branch and re-used (re-positioned) at each step.

function _inversion_parity(seq)
	inv = 0
	for i in 1:length(seq)-1
		for j in i+1:length(seq)
			inv += seq[i] > seq[j]
		end
	end
	return isodd(inv) ? -1 : 1
end

_z2space(d0::Int, d1::Int) = (d1 == 0) ? Z2Space(0=>d0) : (d0 == 0) ? Z2Space(1=>d1) : Z2Space(0=>d0, 1=>d1)

# svd of a (possibly empty) matrix; the number of kept singular values is
# reduced to the numerical rank (σ > tol·σmax)
function _rank_svd(M::AbstractMatrix, tol::Real=1.0e-14)
	T = eltype(M)
	(minimum(size(M)) == 0) && return Matrix{T}(undef, size(M, 1), 0), Float64[], Matrix{T}(undef, size(M, 2), 0)
	U, S, V = LinearAlgebra.svd(M)
	σmax = maximum(S)
	(σmax == 0) && return Matrix{T}(undef, size(M, 1), 0), Float64[], Matrix{T}(undef, size(M, 2), 0)
	keep = S .> σmax * tol
	return U[:, keep], S[keep], V[:, keep]
end

"""
	_fockpropagator_sites(fockstate, brel, krel, L; δ) -> Vector{MPSTensor}

Build the site tensors of the propagator GrassmannMPS on a span of `L`
sites. `brel`/`krel` are the 0-based positions (within the span) of the
bra (conjugated) resp. ket Grassmann variables of the `M = length(brel)`
bands; span sites which carry no variable stay in the vacuum sector.

The coherent-state coefficient tensor is factorized by a left-to-right
SVD sweep over only the `2M` sites carrying a Grassmann variable (cost
O(4^M), independent of the span); the spectator sites in between are
filled with identity pass-through tensors.
"""
# identity pass-through tensor on a spectator site with bond space V:
# only the x=0 sector of the physical index is populated
function _wire_tensor(V::Z2Space, T::Type)
	data = Dict{Z2Irrep, Matrix{T}}()
	for c in (Z2Irrep(0), Z2Irrep(1))
		d0 = dim(V, c)
		d0 == 0 && continue
		d1 = dim(V, Z2Irrep(1 - c.n))
		m = zeros(T, d0 + d1, d0)
		for i in 1:d0
			m[i, i] = one(T)
		end
		data[c] = m
	end
	return TensorMap(data, V ⊗ _ph, V)
end

function _fockpropagator_sites(fockstate::AbstractMatrix, brel::Vector{Int}, krel::Vector{Int},
								L::Int; δ::Float64=1.0e-10)
	M = length(brel)
	(length(krel) == M) || throw(DimensionMismatch("bra and ket positions do not match"))
	allunique(vcat(brel, krel)) || throw(ArgumentError("bra and ket grassmann variables overlap"))
	all(0 .<= vcat(brel, krel) .<= L-1) || throw(ArgumentError("variable outside of the span"))
	T = eltype(fockstate)

	# active sites (carrying a grassmann variable) in span order
	act = sort(vcat(brel, krel))
	Lact = length(act)
	actind = Dict(p => j for (j, p) in enumerate(act))
	# which grassmann variable is carried by each active site
	avars = Vector{Tuple{Symbol, Int}}(undef, Lact)
	for i in 1:M
		avars[actind[brel[i]]] = (:bra, i)
		avars[actind[krel[i]]] = (:ket, i)
	end
	occ(b, i) = (b >> (M - i)) & 1

	# dense coefficient vector over the active sites; the occupation bitstring
	# (y_1, ..., y_{2M}) is encoded as Σ_j y_j·2^(j-1) (y_1 = LSB), matching
	# the Z2 graded basis ordering where the first tensor index varies fastest
	Tvec = zeros(T, 2^Lact)
	for row in 1:2^M, col in 1:2^M
		v = fockstate[row, col]
		abs(v) < δ && continue
		m, n = row - 1, col - 1
		if count(i -> occ(m, i) == 1, 1:M) != count(i -> occ(n, i) == 1, 1:M)
			@warn "Ignore non-physical element: $row $col => $v"
			continue
		end
		ys = zeros(Int, Lact)
		for j in 1:Lact
			ys[j] = avars[j][1] === :bra ? occ(m, avars[j][2]) : occ(n, avars[j][2])
		end
		yint = sum(ys[j] * 2^(j-1) for j in 1:Lact)
		iseven(count_ones(yint)) || throw(ArgumentError("odd number of grassmann variables"))
		# sign from reordering c̄(bra, ascending bands) c(ket, descending bands) into site order
		seq = vcat([brel[i] for i in 1:M if occ(m, i) == 1],
				   [krel[i] for i in M:-1:1 if occ(n, i) == 1])
		Tvec[yint+1] += _inversion_parity(seq) * v
	end
	all(iszero, Tvec) && throw(ArgumentError("empty propagator"))

	# Z2-resolved SVD sweep over the active sites, left to right. The carry R
	# has rows grouped as [even bond basis; odd bond basis] and columns
	# indexed by the integer encoding of the remaining occupations; the next
	# active site is the least significant bit of the column index
	asites = Vector{Any}(undef, Lact)
	d_even, d_odd = 1, 0
	R = reshape(Tvec, 1, 2^Lact)
	for k in 1:Lact-1
		C = size(R, 2)
		# split off y_k (the LSB of the column index): rows become the fused
		# (y_k, bond) basis with y_k major (bond index varies fastest)
		R2 = vcat(R[:, 1:2:C], R[:, 2:2:C])
		d = d_even + d_odd
		# rows of the Z2Irrep(0) block: (y_k=0, l even) followed by (y_k=1, l odd)
		erows = vcat(1:d_even, d+d_even+1:2d)
		# rows of the Z2Irrep(1) block: (y_k=0, l odd) followed by (y_k=1, l even)
		orows = vcat(d_even+1:d, d+1:d+d_even)
		ncol = C ÷ 2
		ecols = [j for j in 1:ncol if iseven(count_ones(j-1))]
		ocols = [j for j in 1:ncol if isodd(count_ones(j-1))]

		U0, S0, V0 = _rank_svd(R2[erows, ecols])
		U1, S1, V1 = _rank_svd(R2[orows, ocols])
		r0, r1 = length(S0), length(S1)

		l_space = _z2space(d_even, d_odd)
		r_space = _z2space(r0, r1)
		data = Dict{Z2Irrep, Matrix{T}}()
		(r0 > 0) && (data[Z2Irrep(0)] = U0)
		(r1 > 0) && (data[Z2Irrep(1)] = U1)
		asites[k] = TensorMap(data, l_space ⊗ _ph, r_space)

		# new carry, with the columns indexed by the remaining occupations
		R = zeros(T, r0+r1, ncol)
		(r0 > 0) && (R[1:r0, ecols] = Diagonal(S0) * V0')
		(r1 > 0) && (R[r0+1:r0+r1, ocols] = Diagonal(S1) * V1')
		d_even, d_odd = r0, r1
	end
	# last active site: (bond ⊗ y ← oneunit)
	d = d_even + d_odd
	data = Dict{Z2Irrep, Matrix{T}}()
	data[Z2Irrep(0)] = reshape(vcat(R[1:d_even, 1], R[d_even+1:d, 2]), d, 1)
	asites[Lact] = TensorMap(data, _z2space(d_even, d_odd) ⊗ _ph, oneunit(_ph))

	# assemble the full span: active tensors at their positions, identity
	# pass-through tensors on the spectator sites in between
	sites = Vector{Any}(undef, L)
	ia = 0
	for p in 0:L-1
		if haskey(actind, p)
			ia += 1
			sites[p+1] = asites[ia]
		else
			sites[p+1] = _wire_tensor(dual(space_r(sites[p])), T)  # p>=1 since pmin is active
		end
	end
	return sites
end

"""
	fockpropagator_gmps(fockstate, lattice, idx; branch, δ, cache) -> GrassmannMPS

Directly convert the Fock-space propagator matrix `fockstate` (i.e. the
operator `exp(coeff*H)` in the occupation number basis) into a
GrassmannMPS on `lattice`. The bra (conjugated) variables sit on time
slice `idx+1` and the ket variables on `idx` (swapped for the backward
branch), exactly as in `sysdynamics_util2`.

If `cache` (a `Dict`) is given, the site tensors of the propagator are
cached by their relative layout and re-used for all time steps.
"""
function fockpropagator_gmps(fockstate::AbstractMatrix, lattice::AbstractGrassmannLattice,
							 idx::Int, branch::Symbol; δ::Float64=1.0e-10,
							 cache::Union{Nothing, Dict}=nothing)
	M = lattice.bands
	(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
	(size(fockstate, 1) == 2^M) || throw(DimensionMismatch("fockstate size does not match bands"))

	ib, ik = branch == :- ? (idx, idx+1) : (idx+1, idx)
	bpos = [index(lattice, ib, conj=true, branch=branch, band=i) for i in 1:M]
	kpos = [index(lattice, ik, conj=false, branch=branch, band=i) for i in 1:M]
	pmin = min(minimum(bpos), minimum(kpos))
	pmax = max(maximum(bpos), maximum(kpos))
	L = pmax - pmin + 1

	brel = bpos .- pmin
	krel = kpos .- pmin
	sites = isnothing(cache) ? nothing : get(cache, (branch, brel, krel), nothing)
	if isnothing(sites)
		sites = _fockpropagator_sites(fockstate, brel, krel, L; δ=δ)
		isnothing(cache) || (cache[(branch, brel, krel)] = sites)
	end

	gmps = vacuumstate(lattice)
	if !(scalartype(sites[1]) <: Real) && (scalartype(gmps[1]) <: Real)
		gmps = complex(gmps)
	end
	for (j, p) in enumerate(pmin:pmax)
		gmps[p] = sites[j]
	end
	unset_svectors!(gmps)
	return gmps
end

function sysdynamics_util2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
								idx::Int=1, branch::Symbol=:+, trunc::TruncationScheme=DefaultKTruncation,
								cache::Union{Nothing, Dict}=nothing)
	H = fockmatrix(model, lattice.bands)
	if branch == :+
		coeff = - im * lattice.δt
	elseif branch == :-
		coeff = im * lattice.δt
	elseif branch == :τ
		coeff = - lattice.δτ
	else
		throw(ArgumentError("branch must be one of :+, :- or :τ"))
	end

	# exact propagator from the spectral decomposition of H
	vals, vecs = eigen(H)
	fockstate = vecs * Diagonal(exp.(coeff .* vals)) * vecs'

	state = fockpropagator_gmps(fockstate, lattice, idx, branch; cache=cache)
	return mult(state, gmps, trunc=trunc)
end

function sysdynamics_forward2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
									trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nt
		gmps = sysdynamics_util2_new(gmps, lattice, model; idx=i, branch=:+, trunc=trunc, cache=cache)
	end
	return gmps
end
function sysdynamics_backward2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
									trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nt
		gmps = sysdynamics_util2_new(gmps, lattice, model; idx=i, branch=:-, trunc=trunc, cache=cache)
	end
	return gmps
end
function sysdynamics_imaginary2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
									trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nτ
		gmps = sysdynamics_util2_new(gmps, lattice, model; idx=i, branch=:τ, trunc=trunc, cache=cache)
	end
	return gmps
end

"""
	sysdynamics2_new(lattice, model; branch, trunc) -> GrassmannMPS

Build the GrassmannMPS of the impurity propagator `K` on `lattice`,
identical to `sysdynamics2` but with the Fock-space propagator converted
directly into coherent-state site tensors (see `fockpropagator_gmps`).
"""
function sysdynamics2_new(lattice::ImagGrassmannLattice, model::AbstractImpurityHamiltonian;
							trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	return sysdynamics_imaginary2_new(gmps, lattice, model; trunc=trunc)
end

function sysdynamics2_new(lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian;
							branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	if isnothing(branch)
		gmps = sysdynamics_forward2_new(gmps, lattice, model; trunc=trunc)
		return sysdynamics_backward2_new(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward2_new(gmps, lattice, model; trunc=trunc) :
								sysdynamics_backward2_new(gmps, lattice, model; trunc=trunc)
	end
end

function sysdynamics2_new(lattice::MixedGrassmannLattice, model::AbstractImpurityHamiltonian;
							branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	if isnothing(branch)
		gmps = sysdynamics_forward2_new(gmps, lattice, model; trunc=trunc)
		gmps = sysdynamics_backward2_new(gmps, lattice, model; trunc=trunc)
		return sysdynamics_imaginary2_new(gmps, lattice, model; trunc=trunc)
	else
		if branch == :+
			return sysdynamics_forward2_new(gmps, lattice, model; trunc=trunc)
		elseif branch == :-
			return sysdynamics_backward2_new(gmps, lattice, model; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return sysdynamics_imaginary2_new(gmps, lattice, model; trunc=trunc)
		end
	end
end
