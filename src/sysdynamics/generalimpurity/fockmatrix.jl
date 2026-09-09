# FockMatrix
# ----------
# A matrix in the Fock (occupation number) basis of the impurity, together
# with the number of bands it belongs to. The matrix has dimension
# 2^bands × 2^bands and is ordered as in `jw_operators` / `fockstate`:
# band 1 is the most significant bit.

"""
	FockMatrix(data, bands)

A matrix in the Fock (occupation number) basis of an impurity with
`bands` fermionic bands. `data` must be a square matrix of size
2^bands × 2^bands, with the occupation of band 1 as the most
significant bit of the basis label (the convention of `jw_operators`
and `fockstate`).
"""
struct FockMatrix{T <: Number}
	data::Matrix{T}
	bands::Int

	function FockMatrix(data::AbstractMatrix{T}, bands::Int) where {T <: Number}
		size(data) == (2^bands, 2^bands) || throw(DimensionMismatch("expected a $(2^bands)×$(2^bands) matrix for bands=$bands, got $(size(data))"))
		return new{T}(convert(Matrix{T}, data), bands)
	end
end

# infer the number of bands from the matrix size (must be a power of two)
function FockMatrix(data::AbstractMatrix)
	n = size(data, 1)
	(size(data, 2) == n && n >= 1 && (n & (n - 1)) == 0) || throw(DimensionMismatch("matrix size $(size(data)) is not 2^bands × 2^bands"))
	return FockMatrix(data, round(Int, log2(n)))
end

Base.size(x::FockMatrix) = size(x.data)
Base.getindex(x::FockMatrix, i::Int, j::Int) = getindex(x.data, i, j)
TK.scalartype(::Type{FockMatrix{T}}) where {T <: Number} = T
TK.scalartype(x::FockMatrix) = scalartype(typeof(x))

# internal helpers: exact propagator / thermal state from the Fock-space
# Hamiltonian matrix of a model with `bands` bands
function _fock_propagator(H::AbstractMatrix, bands::Int, branch::Symbol, dt::Real)
	(branch in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
	coeff = branch == :+ ? - im * dt : branch == :- ? im * dt : - dt
	vals, vecs = eigen(H)
	return FockMatrix(vecs * Diagonal(exp.(coeff .* vals)) * vecs', bands)
end

function _fock_thermalstate(H::AbstractMatrix, bands::Int, β::Real)
	vals, vecs = eigen(H)
	if β == Inf
		rho = vecs[:,1] * vecs[:,1]'
	else
		# shift by E0 to avoid overflow at large β, then normalize by the
		# trace so that tr(rho) == 1
		E0 = minimum(vals)
		weights = exp.(-β .* (vals .- E0))
		weights ./= sum(weights)
		rho = vecs * Diagonal(weights) * vecs'
	end
	return FockMatrix(rho, bands)
end

"""
	fock_propagator(h, branch, dt) -> FockMatrix

The impurity propagator `exp(coeff·Ĥ)` in the Fock basis, where the
coefficient is determined by the Keldysh branch and the time step:
`-im·dt` for `:+`, `im·dt` for `:-` and `-dt` for `:τ`. Interface
function of `ConstImpurityHamiltonian`: the model carries its own
`bands` field (see `num_bands`).
"""
fock_propagator(h::ImpurityHamiltonian, branch::Symbol, dt::Real) =
	_fock_propagator(fockmatrix(h, h.bands), h.bands, branch, dt)

fock_propagator(m::AndersonIM, branch::Symbol, dt::Real) =
	_fock_propagator(fockmatrix(m), num_bands(m), branch, dt)

fock_propagator(m::ToulouseIM, branch::Symbol, dt::Real) =
	_fock_propagator(fockmatrix(m), num_bands(m), branch, dt)

# quench protocol: the imaginary-time branch evolves with the pre-quench h0,
# the real-time branches with the post-quench h1
fock_propagator(h::QuenchImpurityHamiltonian, branch::Symbol, dt::Real) =
	_fock_propagator(fockmatrix(branch == :τ ? h.h0 : h.h1, h.bands), h.bands, branch, dt)

"""
	fock_thermalstate(h, β) -> FockMatrix

The normalized thermal equilibrium state `exp(-βĤ)/tr(exp(-βĤ))` of the
impurity in the Fock basis — the second interface function of
`AbstractImpurityHamiltonian`, needed to set the initial state on
real-time lattices. For `β == Inf` the (normalized) ground state
projector is returned.
"""
fock_thermalstate(h::ImpurityHamiltonian, β::Real) =
	_fock_thermalstate(fockmatrix(h, h.bands), h.bands, β)

fock_thermalstate(m::AndersonIM, β::Real) =
	_fock_thermalstate(fockmatrix(m), num_bands(m), β)

fock_thermalstate(m::ToulouseIM, β::Real) =
	_fock_thermalstate(fockmatrix(m), num_bands(m), β)

# the quench starts from the thermal equilibrium of the pre-quench h0
fock_thermalstate(h::QuenchImpurityHamiltonian, β::Real) =
	fock_thermalstate(ImpurityHamiltonian(h.h0, h.bands), β)

# Fock matrix -> SparseGMPS
# -------------------------
# The coherent-state coefficient tensor
#     T[x] = Σ_{m,n} (-1)^{inv(m,n)} ρ[m,n]
# (x = occupations of the bra/ket Grassmann variables in site order,
# ρ = Fock matrix, inv = inversions needed to sort the monomial
#  c̄(bra, ascending bands) c(ket, descending bands) into site order)
# is assembled directly in the Z2Irrep(0) block of a graded tensor and
# then factorized into site tensors by a Z2-graded SVD sweep. Only the
# sites carrying a Grassmann variable are active; the spectator sites in
# between carry identity pass-through tensors, and sites that reduce to
# the trivial unit tensor are dropped from the SparseGMPS.

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

"""
	_fockpropagator_sites(fockstate, brel, krel, L; δ) -> (sites, trivial)

Build the site tensors of the propagator GrassmannMPS on a span of `L`
sites. `brel`/`krel` are the 0-based positions (within the span) of the
bra (conjugated) resp. ket Grassmann variables of the `M = length(brel)`
bands; span sites which carry no variable stay in the vacuum sector.

The coherent-state coefficient tensor is factorized by a left-to-right
SVD sweep over only the `2M` sites carrying a Grassmann variable (cost
O(4^M), independent of the span); the spectator sites in between are
filled with identity pass-through tensors.

Also returned is a parallel `trivial` flag for every span site, decided
structurally during the sweep (no numerical tolerance): an active site
is trivial iff both of its bonds are oneunit and the vacuum sector
passes through with amplitude exactly one; a spectator site is trivial
iff its bond space is oneunit (the pass-through tensor then reduces to
the unit tensor by construction).
"""
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
	atrivial = Vector{Bool}(undef, Lact)
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
		# the site tensor is the unit tensor iff both bonds are oneunit and
		# the (then 1×1) even block passes the vacuum through with amplitude
		# exactly one (exact comparison: the amplitudes live in the carry)
		atrivial[k] = (d_even == 1 && d_odd == 0 && r0 == 1 && r1 == 0 &&
						U0 == fill(one(T), 1, 1))

		# new carry, with the columns indexed by the remaining occupations
		R = zeros(T, r0+r1, ncol)
		(r0 > 0) && (R[1:r0, ecols] = Diagonal(S0) * V0')
		(r1 > 0) && (R[r0+1:r0+r1, ocols] = Diagonal(S1) * V1')
		d_even, d_odd = r0, r1
	end
	# last active site: (bond ⊗ y ← oneunit); its tensor carries the final
	# amplitudes, so it is the unit tensor iff the left bond is oneunit and
	# the remaining vacuum amplitude is exactly one
	d = d_even + d_odd
	data = Dict{Z2Irrep, Matrix{T}}()
	data[Z2Irrep(0)] = reshape(vcat(R[1:d_even, 1], R[d_even+1:d, 2]), d, 1)
	asites[Lact] = TensorMap(data, _z2space(d_even, d_odd) ⊗ _ph, oneunit(_ph))
	atrivial[Lact] = (d_even == 1 && d_odd == 0 && R[1, 1] == one(T))

	# assemble the full span: active tensors at their positions, identity
	# pass-through tensors on the spectator sites in between
	sites = Vector{Any}(undef, L)
	trivial = Vector{Bool}(undef, L)
	ia = 0
	for p in 0:L-1
		if haskey(actind, p)
			ia += 1
			sites[p+1] = asites[ia]
			trivial[p+1] = atrivial[ia]
		else
			V = dual(space_r(sites[p]))  # p>=1 since pmin is active
			sites[p+1] = _wire_tensor(V, T)
			# a pass-through tensor on a oneunit bond is the unit tensor
			trivial[p+1] = isoneunit(V)
		end
	end
	return sites, trivial
end

"""
	_tosparsegmps(lattice, fm, bpos, kpos; δ) -> SparseGMPS

Convert the Fock matrix `fm` (an operator in the occupation number
basis, e.g. a propagator or a density matrix) into a SparseGMPS on
`lattice`, with the bra (conjugated) variables at the site positions
`bpos` and the ket variables at `kpos`. Sites whose tensors are the
trivial unit tensor — decided structurally during the SVD sweep, see
`_fockpropagator_sites` — are dropped from the SparseGMPS.
"""
function _tosparsegmps(lattice::AbstractGrassmannLattice, fm::FockMatrix, bpos::Vector{Int}, kpos::Vector{Int};
						δ::Float64=1.0e-10)
	M = lattice.bands
	(fm.bands == M) || throw(DimensionMismatch("FockMatrix bands $(fm.bands) do not match lattice bands $M"))
	(length(bpos) == M && length(kpos) == M) || throw(DimensionMismatch("bra/ket positions do not match bands"))

	pmin = min(minimum(bpos), minimum(kpos))
	pmax = max(maximum(bpos), maximum(kpos))
	sites, trivial = _fockpropagator_sites(fm.data, bpos .- pmin, kpos .- pmin, pmax - pmin + 1; δ=δ)

	# keep only the nontrivial site tensors
	data = [sites[j] for j in 1:length(sites) if !trivial[j]]
	pos = [pmin + j - 1 for j in 1:length(sites) if !trivial[j]]
	all(trivial) && return SparseGMPS(typeof(sites[1]))
	return SparseGMPS(convert(Vector{typeof(sites[1])}, data), pos)
end
