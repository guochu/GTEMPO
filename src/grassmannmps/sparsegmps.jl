# SparseGMPS
# ----------
# A GrassmannMPS in which most sites carry the trivial unit tensor of the
# GMPS multiplication (the vacuum scalar: trivial bonds, value 1 in the even
# sector of the physical index). Only the nontrivial site tensors and their
# positions are stored; there is no scaling field (scaling is implicitly 1).
#
# The sites listed in `positions` hold the stored tensors; every other site
# of the underlying lattice implicitly carries the unit tensor. The
# constructor requires that the positions are strictly increasing, that
# the leftmost and rightmost auxiliary spaces are trivial, and that the
# bond structure is consistent with trivial sites filling the gaps.

"""
	SparseGMPS(data, positions)

A GrassmannMPS in which most site tensors are the trivial unit tensor
(the vacuum scalar). Only the nontrivial site tensors (`data`) and their
positions (`positions`) are stored. Sites of the underlying lattice that
are not listed in `positions` implicitly carry the unit tensor of the
GMPS multiplication, i.e. the SparseGMPS represents the same state as
the fully materialized GrassmannMPS.

The constructor checks that
- the positions are sorted in ascending order (unsorted input throws an
  `ArgumentError` instead of being silently re-ordered),
- the leftmost and rightmost auxiliary spaces are trivial, and
- consecutive stored tensors have matching bonds, which are moreover
  trivial whenever trivial sites fill the gap in between.
"""
struct SparseGMPS{A <: MPSTensor}
	data::Vector{A}
	positions::Vector{Int}

	function SparseGMPS(data::Vector{A}, positions::Vector{Int}) where {A <: MPSTensor}
		(length(data) == length(positions)) || throw(DimensionMismatch("data and positions must have the same length"))
		issorted(positions) || throw(ArgumentError("positions must be sorted in ascending order"))
		allunique(positions) || throw(ArgumentError("positions must be unique"))
		if !isempty(data)
			isoneunit(space_l(data[1])) || throw(ArgumentError("the leftmost auxiliary space must be trivial"))
			isoneunit(space_r(data[end])) || throw(ArgumentError("the rightmost auxiliary space must be trivial"))
			for i in 1:length(data)-1
				# consecutive stored tensors must be bond-compatible
				(space_l(data[i+1]) == dual(space_r(data[i]))) ||
					throw(SpaceMismatch("tensors $(positions[i]) and $(positions[i+1]) have incompatible bonds"))
				# a gap between stored positions is filled with trivial sites,
				# which requires the bond across the gap to be trivial
				(positions[i+1] == positions[i] + 1 || isoneunit(space_r(data[i]))) ||
					throw(ArgumentError("nontrivial bond across a gap of trivial sites between positions $(positions[i]) and $(positions[i+1])"))
			end
		end
		return new{A}(data, positions)
	end
end

# convenience method for general vectors, delegating to the checking
# inner constructor after conversion
SparseGMPS(data::AbstractVector{A}, positions::AbstractVector{Int}) where {A <: MPSTensor} =
	SparseGMPS(convert(Vector{A}, data), convert(Vector{Int}, positions))

SparseGMPS(::Type{A}) where {A <: MPSTensor} = SparseGMPS(A[], Int[])
SparseGMPS() = SparseGMPS(mpstensortype(grassmannpspacetype(), Float64))

TK.scalartype(::Type{<:SparseGMPS{A}}) where {A <: MPSTensor} = scalartype(A)
TK.scalartype(x::SparseGMPS) = scalartype(typeof(x))

positions(x::SparseGMPS) = x.positions
Base.length(x::SparseGMPS) = length(x.data)
Base.isempty(x::SparseGMPS) = isempty(x.data)
Base.getindex(x::SparseGMPS, i::Int) = x.data[i]
Base.copy(x::SparseGMPS) = SparseGMPS(copy(x.data), copy(x.positions))

"""
the trivial unit tensor of the GMPS multiplication on one site:
trivial bonds, value 1 in the even sector and 0 in the odd sector
"""
function _trivial_site(::Type{T}) where {T <: Number}
	v = zeros(T, oneunit(_ph) ⊗ _ph ← oneunit(_ph))
	copy!(block(v, Z2Irrep(0)), ones(T, 1, 1))
	return v
end

"""
	SparseGMPS(x::GTerm) -> SparseGMPS

Convert a `GTerm` (a Grassmann monomial in sorted site order, with the
reordering sign absorbed into the coefficient) into a SparseGMPS with
bond dimension 1. The bond spaces carry the parity of the variables
occupied so far (alternating between the even and the odd one-dimensional
sector), so each tensor at a GTerm position forces the physical variable
to be occupied (with the amplitude on the first position), and the
spectator sites in between are pass-through tensors on the vacuum sector.
All sites outside the span are the trivial unit tensor.
"""
function SparseGMPS(x::GTerm{N}) where {N}
	T = scalartype(x)
	A = mpstensortype(grassmannpspacetype(), T)
	pos = collect(positions(x))
	pmin, pmax = pos[1], pos[end]

	data = Vector{A}(undef, pmax - pmin + 1)
	parity = 0  # parity of the bond to the left of the current site
	for p in pmin:pmax
		j = findfirst(==(p), pos)
		if isnothing(j)
			# spectator site: pass the bond through in the vacuum sector
			V = _z2space(1 - parity, parity)
			data[p-pmin+1] = TensorMap(Dict(Z2Irrep(parity) => fill(one(T), 1, 1)), V ⊗ _ph, V)
		else
			# occupied site: the parity flips across the tensor
			Vl = _z2space(1 - parity, parity)
			Vr = _z2space(parity, 1 - parity)
			val = (j == 1) ? x.coeff : one(x.coeff)
			data[p-pmin+1] = TensorMap(Dict(Z2Irrep(1 - parity) => fill(val, 1, 1)), Vl ⊗ _ph, Vr)
			parity = 1 - parity
		end
	end
	return SparseGMPS(data, collect(pmin:pmax))
end

function SparseGMPS(x::GTerm{0})
	T = scalartype(x)
	t = _trivial_site(T)
	copy!(block(t, Z2Irrep(0)), fill(x.coeff, 1, 1))
	return SparseGMPS([t], [1])
end

"""
	togmps(y::SparseGMPS, L::Int; T) -> GrassmannMPS

Materialize a SparseGMPS into a full GrassmannMPS of length `L`: the
stored tensors are placed at their positions and all other sites carry
the trivial unit tensor.
"""
function togmps(y::SparseGMPS, L::Int; T::Type=scalartype(y))
	gmps = GrassmannMPS(promote_type(T, scalartype(y)), L)
	for (t, p) in zip(y.data, y.positions)
		gmps[p] = t
	end
	unset_svectors!(gmps)
	return gmps
end

# mult!
# -----
"""
	mult!(x::GrassmannMPS, y::SparseGMPS; trunc, verbosity) -> GrassmannMPS

In-place multiplication of a GrassmannMPS `x` by a SparseGMPS `y`,
storing the result in `x`. Equivalent to
`mult!(x, togmps(y, length(x)))`, but only the sites covered by the
nontrivial tensors of `y` are recomputed: outside this window the site
tensors of `x` are reused unchanged.

Only the bonds inside the window are controlled by the truncation
scheme; the bonds outside the window are inherited from `x`. The Schmidt
values (svectors) are refreshed inside the window only, so they may
become stale elsewhere (the next `canonicalize!` recomputes them).
"""
function mult!(x::GrassmannMPS, y::SparseGMPS; trunc::TruncationScheme=DefaultTruncation, verbosity::Int=0)
	L = length(x)
	(isempty(y.data) || all(p -> 1 <= p <= L, y.positions)) || throw(BoundsError())
	# the unit element: x * 1 = x
	isempty(y.data) && return x

	p₀, p₁ = y.positions[1], y.positions[end]
	T = promote_type(scalartype(x), scalartype(y))
	vac = _trivial_site(T)

	# --- left-to-right QR sweep over the window [p₀, p₁] ---
	# the carry has the domain (x-bond, y-bond); the head at p₀ is the
	# naive product tensor with its left bonds fused by `leftfuser`
	ky = 1  # index of the next y-tensor
	ytensor(i::Int) = (ky <= length(y.data) && y.positions[ky] == i) ? y.data[ky] : vac

	tmp5 = g_fuse(_mult_site(x[p₀], ytensor(p₀)), 3)
	Xl, Yl = space_l(x[p₀]), space_l(ytensor(p₀))
	leftfuser = GrassmannTensorMap(isomorphism(T, fuse(Xl, Yl), Xl ⊗ Yl))
	@tensor tmp4[1,4;5,6] := leftfuser[1,2,3] * tmp5[2,3,4,5,6]

	A = tensormaptype(spacetype(x), 2, 1, T)
	res = Vector{A}(undef, p₁ - p₀ + 1)
	for i in p₀:p₁-1
		q, r = leftorth!(tmp4, alg = QR())
		res[i-p₀+1] = get_data(q)
		_renormalize!(x, get_data(r), false)
		(ky <= length(y.data) && y.positions[ky] == i) && (ky += 1)
		@tensor tmp1[1,5,4;2] := r[1,2,3] * GrassmannTensorMap(ytensor(i+1))[3,4,5]
		@tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * GrassmannTensorMap(x[i+1])[4,5,6]
		tmp4 = g_fuse(tmp2, 2)
	end
	# --- last site: fuse the right boundary bonds ---
	_t4 = get_data(tmp4)
	rightfuser = GrassmannTensorMap(isomorphism(T, space(_t4,3)' ⊗ space(_t4,4)', fuse(space(_t4,3), space(_t4,4))))
	@tensor tmp[1,2;5] := tmp4[1,2,3,4] * rightfuser[3,4,5]
	res[end] = get_data(tmp)

	# --- assemble the new data: x outside the window, res inside ---
	newdata = Vector{A}(undef, L)
	for i in 1:p₀-1
		newdata[i] = x[i]
	end
	newdata[p₀:p₁] = res
	for i in p₁+1:L
		newdata[i] = x[i]
	end
	x′ = GrassmannMPS(newdata, copy(x.svectors), scaling(x))

	# --- local right-to-left SVD sweep with truncation over [p₀, p₁] ---
	for i in p₁:-1:p₀+1
		u, s, v, err = stable_tsvd(GrassmannTensorMap(x′[i]), (1,), (2, 3), trunc=trunc)
		x′[i] = get_data(permute(v, (1,2), (3,)))
		nr = _renormalize!(x′, get_data(s), false)
		(verbosity > 1) && println("SVD truncerror at bond $(i): ", sqrt(err * err / (nr * nr + err * err)))
		u2 = u * s
		@tensor tmp[-1,-2;-3] := GrassmannTensorMap(x′[i-1])[-1,-2,1] * u2[1,-3]
		x′[i-1] = get_data(tmp)
		x′.s[i] = get_data(s)
	end
	_renormalize!(x′, x′[p₀], false)

	destory_copy!(x.data, x′.data)
	copy!(x.svectors, x′.svectors)
	setscaling!(x, scaling(x′))
	return x
end
mult(x::GrassmannMPS, y::SparseGMPS; kwargs...) = mult!(copy(x), y; kwargs...)
