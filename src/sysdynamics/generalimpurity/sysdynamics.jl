# sysdynamics
# -----------
# Exact impurity propagator: the Fock-space propagator matrix
# (FockMatrix, see fock_propagator) is converted directly into the site
# tensors of a SparseGMPS in the coherent-state (Grassmann)
# representation (see `_tosparsegmps`), which is multiplied into the
# accumulating state with the sparse mult!.
#
# The propagator of every time step of a branch has the same FockMatrix,
# so the SparseGMPS is built once per branch and re-used (re-positioned)
# at each step via the sparse mult!.

# the propagator SparseGMPS of time step `idx` on `branch`; the bra
# (conjugated) variables sit on time slice `idx+1` and the ket variables
# on `idx` (swapped for the backward branch). The cache is keyed by the
# branch and the relative
# layout of the variables, since the propagator is the same at every step.
function _propagator_sparsegmps(lattice::AbstractGrassmannLattice, fm::FockMatrix, idx::Int, branch::Symbol,
								cache::Union{Nothing, Dict}, t::Union{Nothing, Real}=nothing)
	M = lattice.bands
	ib, ik = branch == :- ? (idx, idx+1) : (idx+1, idx)
	bpos = [index(lattice, ib, conj=true, branch=branch, band=i) for i in 1:M]
	kpos = [index(lattice, ik, conj=false, branch=branch, band=i) for i in 1:M]
	pmin = min(minimum(bpos), minimum(kpos))
	brel, krel = bpos .- pmin, kpos .- pmin

	if isnothing(cache)
		return _tosparsegmps(lattice, fm, bpos, kpos)
	end
	# the cache key carries the branch time t: for time-dependent models the
	# propagator differs from step to step
	cached = get(cache, (branch, t, brel, krel), nothing)
	if isnothing(cached)
		sparse = _tosparsegmps(lattice, fm, bpos, kpos)
		cache[(branch, t, brel, krel)] = (sparse.data, sparse.positions .- pmin)
		return sparse
	end
	data, relpos = cached
	return SparseGMPS(data, relpos .+ pmin)
end

function _sysdynamics_util(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
								idx::Int=1, branch::Symbol=:+, trunc::TruncationScheme=DefaultKTruncation,
								cache::Union{Nothing, Dict}=nothing)
	dt = branch == :τ ? lattice.δτ : lattice.δt
	# left endpoint of the physical time interval of this step (only used by
	# time-dependent models): the forward branch runs from t = 0 to t = T,
	# the backward branch runs back from t = T to t = 0
	t = (branch == :+) ? (idx - 1) * dt : (branch == :- ? (lattice.Nt - idx) * dt : nothing)
	fm = (branch == :τ) ? fock_propagator(model, branch, dt) : fock_propagator(model, branch, dt, t)
	sparse = _propagator_sparsegmps(lattice, fm, idx, branch, cache, t)
	return mult!(gmps, sparse, trunc=trunc)
end

function sysdynamics_forward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
									trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nt
		gmps = _sysdynamics_util(gmps, lattice, model; idx=i, branch=:+, trunc=trunc, cache=cache)
	end
	return gmps
end
function sysdynamics_backward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
									trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nt
		gmps = _sysdynamics_util(gmps, lattice, model; idx=i, branch=:-, trunc=trunc, cache=cache)
	end
	return gmps
end
function sysdynamics_imaginary!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
									trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nτ
		gmps = _sysdynamics_util(gmps, lattice, model; idx=i, branch=:τ, trunc=trunc, cache=cache)
	end
	return gmps
end

"""
	sysdynamics(lattice, model; branch, trunc) -> GrassmannMPS

Build the GrassmannMPS of the impurity propagator `K` on `lattice`:
the Fock-space propagator (`fock_propagator`) is converted directly
into a SparseGMPS (`_tosparsegmps`) which is multiplied into the
accumulating state with the sparse `mult!`.
"""
function sysdynamics(lattice::ImagGrassmannLattice, model::AbstractImpurityHamiltonian;
							trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
end

function sysdynamics(lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian;
							branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	if isnothing(branch)
		gmps = sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? sysdynamics_forward!(gmps, lattice, model; trunc=trunc) :
								sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	end
end

function sysdynamics(lattice::MixedGrassmannLattice, model::AbstractImpurityHamiltonian;
							branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	if isnothing(branch)
		gmps = sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		gmps = sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
		return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
	else
		if branch == :+
			return sysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		elseif branch == :-
			return sysdynamics_backward!(gmps, lattice, model; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return sysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
		end
	end
end

# the bare propagator: the coherent-state overlap ∏_b (1 + c̄_b c_b)
# between the bra and ket variables of the propagator is removed by
# multiplying with its inverse ∏_b (1 - c̄_b c_b) — the overlap is what
# `bulkconnection` re-applies, so bulkconnection on the bare propagator
# recovers the full propagator exactly. The removal follows the original
# implementation: the propagator is materialized on its (small) window,
# the overlap-removal GTerms are applied, and the result is canonicalized.
# The window content depends only on the relative layout of the variables,
# so it is computed once per branch and re-positioned at each time step.
function _bare_window_sparsegmps(lattice::AbstractGrassmannLattice, fm::FockMatrix,
									bwin::Vector{Int}, kwin::Vector{Int}, Lw::Int)
	# bwin/kwin: 1-based positions of the bra/ket variables inside the window
	sites, _ = _fockpropagator_sites(fm.data, bwin .- 1, kwin .- 1, Lw)
	state = GrassmannMPS(convert(Vector{typeof(sites[1])}, sites))

	# multiply by ∏_b (1 - c̄_b c_b): each GTerm connects the bra variable
	# of band b with the ket variable of band b
	for band in 1:length(bwin)
		apply!(exp(GTerm(bwin[band], kwin[band], coeff=-1)), state)
	end
	canonicalize!(state, alg=Orthogonalize(SVD(), trunc=NoTruncation(), normalize=false))

	# fold the scaling into the site tensors (SparseGMPS carries no scaling
	# field): the window GrassmannMPS of length Lw satisfies
	# (∏T_i)·scaling^Lw = Ô_bare, so multiplying each tensor by scaling
	# gives ∏(T_i·scaling) = Ô_bare exactly — the same operator content in
	# every ordering, with the scaling distributed evenly over the window
	data = copy(state.data)
	s = scaling(state)
	for i in 1:Lw
		data[i] = data[i] * s
	end
	return SparseGMPS(data, collect(1:Lw))
end

function _bare_propagator_sparsegmps(lattice::AbstractGrassmannLattice, fm::FockMatrix, idx::Int,
									branch::Symbol, cache::Union{Nothing, Dict}, t::Union{Nothing, Real}=nothing)
	M = lattice.bands
	ib, ik = branch == :- ? (idx, idx+1) : (idx+1, idx)
	bpos = [index(lattice, ib, conj=true, branch=branch, band=i) for i in 1:M]
	kpos = [index(lattice, ik, conj=false, branch=branch, band=i) for i in 1:M]
	pmin = min(minimum(bpos), minimum(kpos))
	pmax = max(maximum(bpos), maximum(kpos))
	bwin, kwin = bpos .- pmin .+ 1, kpos .- pmin .+ 1

	if isnothing(cache)
		sparse = _bare_window_sparsegmps(lattice, fm, bwin, kwin, pmax - pmin + 1)
		return SparseGMPS(sparse.data, sparse.positions .+ (pmin - 1))
	end
	# the cache key carries the branch time t: for time-dependent models the
	# propagator differs from step to step
	cached = get(cache, (branch, t, bwin, kwin), nothing)
	if isnothing(cached)
		sparse = _bare_window_sparsegmps(lattice, fm, bwin, kwin, pmax - pmin + 1)
		cache[(branch, t, bwin, kwin)] = (sparse.data, sparse.positions)
		return SparseGMPS(sparse.data, sparse.positions .+ (pmin - 1))
	end
	data, relpos = cached
	return SparseGMPS(data, relpos .+ (pmin - 1))
end

function _baresysdynamics_util(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
									idx::Int=1, branch::Symbol=:+, trunc::TruncationScheme=DefaultKTruncation,
									cache::Union{Nothing, Dict}=nothing)
	dt = branch == :τ ? lattice.δτ : lattice.δt
	# left endpoint of the physical time interval of this step (only used by
	# time-dependent models): the forward branch runs from t = 0 to t = T,
	# the backward branch runs back from t = T to t = 0
	t = (branch == :+) ? (idx - 1) * dt : (branch == :- ? (lattice.Nt - idx) * dt : nothing)
	fm = (branch == :τ) ? fock_propagator(model, branch, dt) : fock_propagator(model, branch, dt, t)
	sparse = _bare_propagator_sparsegmps(lattice, fm, idx, branch, cache, t)
	return mult!(gmps, sparse, trunc=trunc)
end

"""
	baresysdynamics_forward!(gmps, lattice, model; trunc) -> GrassmannMPS

In-place: multiply the bare propagators (without the coherent-state
bra-ket overlaps) of all time steps of the forward branch into `gmps`.
"""
function baresysdynamics_forward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
										trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nt
		gmps = _baresysdynamics_util(gmps, lattice, model; idx=i, branch=:+, trunc=trunc, cache=cache)
	end
	return gmps
end
"""
	baresysdynamics_backward!(gmps, lattice, model; trunc) -> GrassmannMPS

In-place: multiply the bare propagators of all time steps of the
backward branch into `gmps`.
"""
function baresysdynamics_backward!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
										trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nt
		gmps = _baresysdynamics_util(gmps, lattice, model; idx=i, branch=:-, trunc=trunc, cache=cache)
	end
	return gmps
end
"""
	baresysdynamics_imaginary!(gmps, lattice, model; trunc) -> GrassmannMPS

In-place: multiply the bare propagators of all imaginary-time steps
into `gmps`.
"""
function baresysdynamics_imaginary!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
										trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nτ
		gmps = _baresysdynamics_util(gmps, lattice, model; idx=i, branch=:τ, trunc=trunc, cache=cache)
	end
	return gmps
end

"""
	baresysdynamics(lattice, model; branch, trunc) -> GrassmannMPS

Exact version of `baresysdynamics`: the impurity propagator without the
coherent-state bra-ket overlaps (see `_bare_window_sparsegmps`).
Applying `bulkconnection` to its output gives `sysdynamics`, exactly
as `bulkconnection` on `baresysdynamics` gives `sysdynamics`.
"""
function baresysdynamics(lattice::ImagGrassmannLattice, model::AbstractImpurityHamiltonian;
								trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	return baresysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
end

function baresysdynamics(lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian;
								branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	if isnothing(branch)
		gmps = baresysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		return baresysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? baresysdynamics_forward!(gmps, lattice, model; trunc=trunc) :
								baresysdynamics_backward!(gmps, lattice, model; trunc=trunc)
	end
end

function baresysdynamics(lattice::MixedGrassmannLattice, model::AbstractImpurityHamiltonian;
								branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	if isnothing(branch)
		gmps = baresysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		gmps = baresysdynamics_backward!(gmps, lattice, model; trunc=trunc)
		return baresysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
	else
		if branch == :+
			return baresysdynamics_forward!(gmps, lattice, model; trunc=trunc)
		elseif branch == :-
			return baresysdynamics_backward!(gmps, lattice, model; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return baresysdynamics_imaginary!(gmps, lattice, model; trunc=trunc)
		end
	end
end

# sysdynamics_fast
# ---------------------
# Fast version of sysdynamics. The single-step propagator is a
# SparseGMPS supported on a small window of sites. Whenever the windows
# of all time steps are pairwise disjoint in the ordering of the lattice
# — true for the time-local orderings such as A1B1B̄1Ā1 — the product of
# all propagators is obtained by filling the window tensors directly
# into a vacuum GrassmannMPS: no GMPS multiplications at all, and the
# result is the exact product (no truncation). For general orderings the
# windows overlap; the propagator is then built in a canonical ordering
# in which the windows are disjoint and transformed to the requested
# ordering with changeordering (or, when only windows of *different*
# branches overlap, per-branch GMPSs are tiled and multiplied).

# the canonical ordering in which the propagator windows of one branch
# (for real time: of both branches together) are pairwise disjoint
_fastordering(lattice::ImagGrassmannLattice) = A1B1B̄1Ā1()
_fastordering(lattice::RealGrassmannLattice) = LayoutStyle(lattice) isa TimeLocalLayout ?
		A1B1ā1b̄1Ā1B̄1a1b1() : A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2()
_fastordering(lattice::MixedGrassmannLattice) = A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1()

# the propagator windows of all time steps of the given branches, in the
# ordering of `lattice`; the per-branch window layout is cached, so each
# distinct window is built only once
function _collect_windows(lattice::AbstractGrassmannLattice, model, branches; bare::Bool=false)
	windows = SparseGMPS[]
	for branch in branches
		N = branch == :τ ? lattice.Nτ : lattice.Nt
		N == 0 && continue
		dt = branch == :τ ? lattice.δτ : lattice.δt
		fm = fock_propagator(model, branch, dt)
		cache = Dict{Tuple, Any}()
		for i in 1:N
			sparse = bare ? _bare_propagator_sparsegmps(lattice, fm, i, branch, cache) :
							_propagator_sparsegmps(lattice, fm, i, branch, cache)
			push!(windows, sparse)
		end
	end
	return windows
end

"""
	_try_tile(lattice, windows) -> Union{GrassmannMPS, Nothing}

Fill the window SparseGMPSs directly into a vacuum GrassmannMPS. The
result represents the product of all window operators exactly, provided
their nontrivial positions are pairwise disjoint (the site-wise product
with the implicit vacuum unit tensors then equals the operator product);
otherwise `nothing` is returned.
"""
function _try_tile(lattice::AbstractGrassmannLattice, windows::Vector{<:SparseGMPS})
	allpos = Int[]
	for w in windows
		append!(allpos, w.positions)
	end
	allunique(allpos) || return nothing
	(isempty(allpos) || (minimum(allpos) >= 1 && maximum(allpos) <= length(lattice))) || throw(BoundsError())

	T = scalartype(lattice)
	for w in windows
		T = promote_type(T, scalartype(w))
	end
	gmps = (T == scalartype(lattice)) ? vacuumstate(lattice) : GrassmannMPS(T, length(lattice))
	for w in windows
		for (t, p) in zip(w.data, w.positions)
			gmps[p] = (scalartype(t) == T) ? t : complex(t)
		end
	end
	unset_svectors!(gmps)
	return gmps
end

function _fast_driver(lattice::AbstractGrassmannLattice{O}, model, branches;
							bare::Bool=false, trunc::TruncationScheme=DefaultKTruncation) where {O}
	branches = Tuple(b for b in branches if (b == :τ ? lattice.Nτ : lattice.Nt) > 0)
	isempty(branches) && return vacuumstate(lattice)

	# 1) windows of all branches pairwise disjoint in the target ordering:
	#    direct tiling, exact product of all propagators
	windows = _collect_windows(lattice, model, branches; bare=bare)
	gmps = _try_tile(lattice, windows)
	isnothing(gmps) || return gmps

	# 2) windows of each single branch disjoint (only windows of different
	#    branches overlap): tile per branch and multiply
	if length(branches) > 1
		gs = [_try_tile(lattice, _collect_windows(lattice, model, (b,); bare=bare)) for b in branches]
		if !any(isnothing, gs)
			gmps = gs[1]
			for i in 2:length(gs)
				gmps = mult(gmps, gs[i], trunc=trunc)
			end
			return gmps
		end
	end

	# 3) build in the canonical ordering (windows disjoint) and change
	#    the ordering of the result
	lattice2 = similar(lattice, ordering=_fastordering(lattice))
	windows2 = _collect_windows(lattice2, model, branches; bare=bare)
	gmps2 = _try_tile(lattice2, windows2)
	if isnothing(gmps2)
		# windows of different branches overlap even in the canonical
		# ordering (mixed contour): tile per branch and multiply
		gs2 = [_try_tile(lattice2, _collect_windows(lattice2, model, (b,); bare=bare)) for b in branches]
		any(isnothing, gs2) && error("propagator windows overlap even in the canonical ordering")
		gmps2 = gs2[1]
		for i in 2:length(gs2)
			gmps2 = mult(gmps2, gs2[i], trunc=trunc)
		end
	end
	return changeordering(O, lattice2, gmps2, trunc=trunc)[2]
end

"""
	sysdynamics_fast(lattice, model; branch, trunc) -> GrassmannMPS

Fast version of `sysdynamics`: if the single-step propagator
windows are pairwise disjoint in the ordering of `lattice` (time-local
orderings such as A1B1B̄1Ā1), their tensors are filled directly into a
vacuum GrassmannMPS — the exact product of all propagators, without any
GMPS multiplication. Otherwise the propagator is built in a canonical
ordering (A1B1B̄1Ā1 for imaginary time, A1B1ā1b̄1Ā1B̄1a1b1 resp.
A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2 for real time,
A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1 for mixed contours) and transformed to
the requested ordering with `changeordering`.
"""
function sysdynamics_fast(lattice::ImagGrassmannLattice, model::ConstImpurityHamiltonian;
								trunc::TruncationScheme=DefaultKTruncation)
	return _fast_driver(lattice, model, (:τ,); trunc=trunc)
end

function sysdynamics_fast(lattice::RealGrassmannLattice, model::ConstImpurityHamiltonian;
								branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		branches = (:+, :-)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		branches = (branch,)
	end
	return _fast_driver(lattice, model, branches; trunc=trunc)
end

function sysdynamics_fast(lattice::MixedGrassmannLattice, model::ConstImpurityHamiltonian;
								branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		branches = (:+, :-, :τ)
	elseif branch in (:+, :-, :τ)
		branches = (branch,)
	else
		throw(ArgumentError("branch must be one of :+, :- or :τ"))
	end
	return _fast_driver(lattice, model, branches; trunc=trunc)
end

"""
	baresysdynamics_fast(lattice, model; branch, trunc) -> GrassmannMPS

Fast version of `baresysdynamics`, built by direct tiling of the
bare propagator windows (see `sysdynamics_fast`). Applying
`bulkconnection` to its output gives `sysdynamics_fast`.
"""
function baresysdynamics_fast(lattice::ImagGrassmannLattice, model::ConstImpurityHamiltonian;
									trunc::TruncationScheme=DefaultKTruncation)
	return _fast_driver(lattice, model, (:τ,); bare=true, trunc=trunc)
end

function baresysdynamics_fast(lattice::RealGrassmannLattice, model::ConstImpurityHamiltonian;
									branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		branches = (:+, :-)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		branches = (branch,)
	end
	return _fast_driver(lattice, model, branches; bare=true, trunc=trunc)
end

function baresysdynamics_fast(lattice::MixedGrassmannLattice, model::ConstImpurityHamiltonian;
									branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	if isnothing(branch)
		branches = (:+, :-, :τ)
	elseif branch in (:+, :-, :τ)
		branches = (branch,)
	else
		throw(ArgumentError("branch must be one of :+, :- or :τ"))
	end
	return _fast_driver(lattice, model, branches; bare=true, trunc=trunc)
end
