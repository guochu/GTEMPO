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

# quench protocol: the imaginary-time branch evolves with the pre-quench hτ,
# the real-time branches with the post-quench ht
fock_propagator(h::QuenchedImpurityHamiltonian, branch::Symbol, dt::Real) =
	_fock_propagator(fockmatrix(branch == :τ ? h.hτ : h.ht, h.bands), h.bands, branch, dt)

# time-dependent model: the τ branch evolves with the constant hτ
function fock_propagator(h::TdImpurityHamiltonian, branch::Symbol, dt::Real)
	(branch == :τ) || throw(ArgumentError("the propagator of a time-dependent model on the real-time branches requires the branch time t: use fock_propagator(model, branch, dt, t)"))
	return _fock_propagator(fockmatrix(h.hτ, h.bands), h.bands, :τ, dt)
end

# time-dependent model: the real-time branches evolve with ht + Σ htt(t)
function fock_propagator(h::TdImpurityHamiltonian, branch::Symbol, dt::Real, t::Real)
	(branch in (:+, :-)) || throw(ArgumentError("branch must be one of :+ or :- for the time-dependent propagator"))
	return _fock_propagator(fockmatrix(h(t), h.bands), h.bands, branch, dt)
end

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

# the quench starts from the thermal equilibrium of the pre-quench hτ
fock_thermalstate(h::QuenchedImpurityHamiltonian, β::Real) =
	_fock_thermalstate(fockmatrix(h.hτ, h.bands), h.bands, β)

# the thermal state of a time-dependent model is that of the constant hτ
fock_thermalstate(h::TdImpurityHamiltonian, β::Real) =
	_fock_thermalstate(fockmatrix(h.hτ, h.bands), h.bands, β)

# Fock matrix -> SparseGMPS

# Fock matrix -> SparseGMPS 旧装配路径（_tosparsegmps/_wire_tensor 等）已随
# 旧栈退役（D2 裁决：K 构造统一走 cell automaton，见 src/port/system.jl）。
