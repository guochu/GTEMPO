# sysdynamics2_new
# ---------------
# Same propagator as sysdynamics2, but built from the Fock-space
# propagator matrix (FockMatrix) which is converted directly into the
# site tensors of a SparseGMPS in the coherent-state (Grassmann)
# representation (see `_tosparsegmps`), instead of being decomposed into
# GTerms which are then applied one by one to the vacuum state.
#
# The propagator of every time step of a branch has the same FockMatrix,
# so the SparseGMPS is built once per branch and re-used (re-positioned)
# at each step via the sparse mult!.

# the propagator SparseGMPS of time step `idx` on `branch`; the bra
# (conjugated) variables sit on time slice `idx+1` and the ket variables
# on `idx` (swapped for the backward branch), exactly as in
# `sysdynamics_util2`. The cache is keyed by the branch and the relative
# layout of the variables, since the propagator is the same at every step.
function _propagator_sparsegmps(lattice::AbstractGrassmannLattice, fm::FockMatrix, idx::Int, branch::Symbol,
								cache::Union{Nothing, Dict})
	M = lattice.bands
	ib, ik = branch == :- ? (idx, idx+1) : (idx+1, idx)
	bpos = [index(lattice, ib, conj=true, branch=branch, band=i) for i in 1:M]
	kpos = [index(lattice, ik, conj=false, branch=branch, band=i) for i in 1:M]
	pmin = min(minimum(bpos), minimum(kpos))
	brel, krel = bpos .- pmin, kpos .- pmin

	if isnothing(cache)
		return _tosparsegmps(lattice, fm, bpos, kpos)
	end
	cached = get(cache, (branch, brel, krel), nothing)
	if isnothing(cached)
		sparse = _tosparsegmps(lattice, fm, bpos, kpos)
		cache[(branch, brel, krel)] = (sparse.data, sparse.positions .- pmin)
		return sparse
	end
	data, relpos = cached
	return SparseGMPS(data, relpos .+ pmin)
end

function sysdynamics_util2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
								idx::Int=1, branch::Symbol=:+, trunc::TruncationScheme=DefaultKTruncation,
								cache::Union{Nothing, Dict}=nothing)
	dt = branch == :τ ? lattice.δτ : lattice.δt
	fm = fock_propagator(model, branch, dt, lattice.bands)
	sparse = _propagator_sparsegmps(lattice, fm, idx, branch, cache)
	return mult!(gmps, sparse, trunc=trunc)
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
identical to `sysdynamics2` but with the Fock-space propagator
(`fock_propagator`) converted directly into a SparseGMPS
(`_tosparsegmps`) which is multiplied into the accumulating state with
the sparse `mult!`.
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
									bwin::Vector{Int}, kwin::Vector{Int}, Lw::Int, L::Int)
	# bwin/kwin: 1-based positions of the bra/ket variables inside the
	# window; L: total lattice length
	sites, _ = _fockpropagator_sites(fm.data, bwin .- 1, kwin .- 1, Lw)
	state = GrassmannMPS(convert(Vector{typeof(sites[1])}, sites))

	# multiply by ∏_b (1 - c̄_b c_b): each GTerm connects the bra variable
	# of band b with the ket variable of band b
	for band in 1:length(bwin)
		apply!(exp(GTerm(bwin[band], kwin[band], coeff=-1)), state)
	end
	canonicalize!(state, alg=Orthogonalize(SVD(), trunc=NoTruncation(), normalize=false))

	# fold the scaling into the site tensors (SparseGMPS carries no scaling
	# field). The scaling of the window GrassmannMPS of length Lw satisifies
	# norm = sqrt(⟨ψ|ψ⟩)·scaling^Lw, while the mult! of the full lattice
	# expects norm = sqrt(⟨ψ|ψ⟩)·(per-site factor)^L, so each window tensor
	# carries scaling^(L/Lw) — the same total as scaling^L, distributed as
	# in GrassmannMPS (which stores scaling^(1/L) per site via _rescaling!)
	data = copy(state.data)
	s = scaling(state)^(L / Lw)
	for i in 1:Lw
		data[i] = data[i] * s
	end
	return SparseGMPS(data, collect(1:Lw))
end

function _bare_propagator_sparsegmps(lattice::AbstractGrassmannLattice, fm::FockMatrix, idx::Int,
									branch::Symbol, cache::Union{Nothing, Dict})
	M = lattice.bands
	ib, ik = branch == :- ? (idx, idx+1) : (idx+1, idx)
	bpos = [index(lattice, ib, conj=true, branch=branch, band=i) for i in 1:M]
	kpos = [index(lattice, ik, conj=false, branch=branch, band=i) for i in 1:M]
	pmin = min(minimum(bpos), minimum(kpos))
	pmax = max(maximum(bpos), maximum(kpos))
	bwin, kwin = bpos .- pmin .+ 1, kpos .- pmin .+ 1

	if isnothing(cache)
		sparse = _bare_window_sparsegmps(lattice, fm, bwin, kwin, pmax - pmin + 1, length(lattice))
		return SparseGMPS(sparse.data, sparse.positions .+ (pmin - 1))
	end
	cached = get(cache, (branch, bwin, kwin), nothing)
	if isnothing(cached)
		sparse = _bare_window_sparsegmps(lattice, fm, bwin, kwin, pmax - pmin + 1, length(lattice))
		cache[(branch, bwin, kwin)] = (sparse.data, sparse.positions)
		return SparseGMPS(sparse.data, sparse.positions .+ (pmin - 1))
	end
	data, relpos = cached
	return SparseGMPS(data, relpos .+ (pmin - 1))
end

function baresysdynamics_util2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
									idx::Int=1, branch::Symbol=:+, trunc::TruncationScheme=DefaultKTruncation,
									cache::Union{Nothing, Dict}=nothing)
	dt = branch == :τ ? lattice.δτ : lattice.δt
	fm = fock_propagator(model, branch, dt, lattice.bands)
	sparse = _bare_propagator_sparsegmps(lattice, fm, idx, branch, cache)
	return mult!(gmps, sparse, trunc=trunc)
end

function baresysdynamics_forward2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
										trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nt
		gmps = baresysdynamics_util2_new(gmps, lattice, model; idx=i, branch=:+, trunc=trunc, cache=cache)
	end
	return gmps
end
function baresysdynamics_backward2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
										trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nt
		gmps = baresysdynamics_util2_new(gmps, lattice, model; idx=i, branch=:-, trunc=trunc, cache=cache)
	end
	return gmps
end
function baresysdynamics_imaginary2_new(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, model;
										trunc::TruncationScheme=DefaultKTruncation)
	cache = Dict{Tuple, Any}()
	for i in 1:lattice.Nτ
		gmps = baresysdynamics_util2_new(gmps, lattice, model; idx=i, branch=:τ, trunc=trunc, cache=cache)
	end
	return gmps
end

"""
	baresysdynamics2_new(lattice, model; branch, trunc) -> GrassmannMPS

Exact version of `baresysdynamics`: the impurity propagator without the
coherent-state bra-ket overlaps (see `_bare_window_sparsegmps`).
Applying `bulkconnection` to its output gives `sysdynamics2_new`, exactly
as `bulkconnection` on `baresysdynamics` gives `sysdynamics`.
"""
function baresysdynamics2_new(lattice::ImagGrassmannLattice, model::AbstractImpurityHamiltonian;
								trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	return baresysdynamics_imaginary2_new(gmps, lattice, model; trunc=trunc)
end

function baresysdynamics2_new(lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian;
								branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	if isnothing(branch)
		gmps = baresysdynamics_forward2_new(gmps, lattice, model; trunc=trunc)
		return baresysdynamics_backward2_new(gmps, lattice, model; trunc=trunc)
	else
		(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :-"))
		return (branch == :+) ? baresysdynamics_forward2_new(gmps, lattice, model; trunc=trunc) :
								baresysdynamics_backward2_new(gmps, lattice, model; trunc=trunc)
	end
end

function baresysdynamics2_new(lattice::MixedGrassmannLattice, model::AbstractImpurityHamiltonian;
								branch::Union{Nothing, Symbol}=nothing, trunc::TruncationScheme=DefaultKTruncation)
	gmps = vacuumstate(lattice)
	if isnothing(branch)
		gmps = baresysdynamics_forward2_new(gmps, lattice, model; trunc=trunc)
		gmps = baresysdynamics_backward2_new(gmps, lattice, model; trunc=trunc)
		return baresysdynamics_imaginary2_new(gmps, lattice, model; trunc=trunc)
	else
		if branch == :+
			return baresysdynamics_forward2_new(gmps, lattice, model; trunc=trunc)
		elseif branch == :-
			return baresysdynamics_backward2_new(gmps, lattice, model; trunc=trunc)
		else
			(branch == :τ) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
			return baresysdynamics_imaginary2_new(gmps, lattice, model; trunc=trunc)
		end
	end
end
