# shared utilities for all test parts: spectra, Grassmann orderings and
# exact-diagonalization (ED) reference builders in the fermionic bit basis

# ---------------------------------------------------------------- spectra
f(D, ϵ) = sqrt(1-(ϵ/D)^2) / π
spectrum_func(D) = spectrum(ϵ->f(D, ϵ), lb=-D, ub=D)   # semicircular spectrum
spectrum_func() = spectrum_func(10)
spectrum_func2(D=10) = spectrum(ϵ->f(D, ϵ), lb=0, ub=D) # half-line spectrum

# ------------------------------------------------------- Grassmann orderings
const imag_orderings   = [A1Ā1B1B̄1(), A1B1B̄1Ā1()]
const real_orderings   = [A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(),
						  A1B1ā1b̄1Ā1B̄1a1b1(),
						  A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(),
						  A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2()]
const mixed_orderings  = [A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2(),
						  A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2(),
						  A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2()]

# ---------------------------------------------------------------- metrics
relerr(a, b) = norm(a - b) / norm(a)

function _error(a::Number, b::Number, tol::Real)
	ab = a - b
	abs(a) < tol ? abs(ab) : abs(ab / a)
end

# ------------------------------------------------- single-particle matrices
using LinearAlgebra: Diagonal, eigen, Hermitian, tr

function spin_half_matrices()
	s_SP = Array{Float64, 2}([0 0; 1 0])
	s_SM = Array{Float64, 2}([0 1; 0 0])
	s_Z = Array{Float64, 2}([-1 0; 0 1])
	s_x = s_SP+s_SM
	s_y = -im*(s_SP-s_SM)
	n = Array{Float64, 2}([0 0; 0 1])
	return Dict("x"=>s_x, "y"=>s_y, "z"=>s_Z, "+"=>s_SP, "-"=>s_SM, "n"=>n)
end

function Aop(d::Int)
	(d <= 1) && error("d must be larger than 1.")
	a = zeros(Float64, d, d)
	for i = 1:(d - 1)
		a[i, i+1] = sqrt(i)
	end
	return a
end

ADAGop(d::Int) = Array(transpose(Aop(d)))
Nop(d::Int) = ADAGop(d) * Aop(d)

function boson(;d::Int=5)
	_N = Nop(d)
	_N2 = _N * _N
	return Dict("a"=>Aop(d),"adag"=>ADAGop(d), "n"=>_N, "n2"=>_N2)
end

# --------------------------------------- ED in the fermionic bit basis
# mode m (1-based) corresponds to bit (m-1); Jordan-Wigner sign
_fermion_sign(state::Int, m::Int) =
	iseven(count_ones(state & ((1 << (m-1)) - 1))) ? 1 : -1

# apply a sequence of creation/annihilation ops (applied left to right:
# ops[1] acts first) on the Fock state `state` (bitmask); returns (sign, state)
# with sign=0 if the process is forbidden
function _apply_fermion_seq(state::Int, ops::Vector{Tuple{Int,Bool}})
	sign = 1
	s = state
	for (m, create) in ops
		bit = 1 << (m-1)
		occupied = (s & bit) != 0
		if create
			occupied && return 0, 0
			sign *= _fermion_sign(s, m)
			s |= bit
		else
			occupied || return 0, 0
			sign *= _fermion_sign(s, m)
			s &= ~bit
		end
	end
	return sign, s
end

# add coeff * c†_i c_j  (tunneling)
function ed_tunneling!(H::AbstractMatrix, i::Int, j::Int, coeff=1.0)
	dim = size(H, 1)
	for s in 0:dim-1
		sign, t = _apply_fermion_seq(s, [(j, false), (i, true)])
		(sign != 0) && (H[t+1, s+1] += coeff * sign)
	end
	return H
end

# add coeff * n_i  (density)
ed_density!(H::AbstractMatrix, i::Int, coeff=1.0) = ed_tunneling!(H, i, i, coeff)

# add coeff * c†_i c†_j c_k c_l  (interaction; ops applied as a_l first)
function ed_interaction!(H::AbstractMatrix, i::Int, j::Int, k::Int, l::Int, coeff=1.0)
	dim = size(H, 1)
	for s in 0:dim-1
		sign, t = _apply_fermion_seq(s, [(l, false), (k, false), (j, true), (i, true)])
		(sign != 0) && (H[t+1, s+1] += coeff * sign)
	end
	return H
end

# add coeff * n_i n_j  (density-density interaction)
function ed_density_interaction!(H::AbstractMatrix, i::Int, j::Int, coeff=1.0)
	dim = size(H, 1)
	for s in 0:dim-1
		occ_i, occ_j = ((s >> (i-1)) & 1) == 1, ((s >> (j-1)) & 1) == 1
		(occ_i && occ_j) && (H[s+1, s+1] += coeff)
	end
	return H
end

# annihilators c[1..m] as dense matrices in the 2^m-dimensional bit basis
function fermion_annihilators(m::Int)
	dim = 1 << m
	cs = [zeros(Float64, dim, dim) for _ in 1:m]
	for i in 1:m, s in 0:dim-1
		sign, t = _apply_fermion_seq(s, [(i, false)])
		(sign != 0) && (cs[i][t+1, s+1] = sign)
	end
	return cs
end

# ------------------------------------------------------- ED model builders
# Fermionic impurity coupled to fermionic single-mode bath modes:
#   H = μ Σ_imp n_i + U n₁n₂ + Σ_m [ω_m b†_m b_m + √α_m (a†_m b_m + b†_m a_m)]
# `bathspecs` lists one (ω, α) per bath mode; `nbands` is the number of
# impurity bands (default: one bath mode per band). Impurity modes are
# 1..nbands, bath modes nbands+1..nbands+length(bathspecs).
# Returns (H, a, adag, H0) where H0 is the decoupled Hamiltonian (impurity +
# bath only) used to build separable thermal initial states, and a/adag
# annihilate on impurity band 1.
function singlemode_ed(; μ::Real, U::Real, bathspecs::Vector{<:Tuple{Real,Real}}, nbands::Int=length(bathspecs))
	nmodes = nbands + length(bathspecs)
	dim = 1 << nmodes
	H = zeros(Float64, dim, dim)
	for i in 1:nbands
		ed_density!(H, i, μ)
	end
	if (nbands == 2) && (U != zero(U))
		ed_density_interaction!(H, 1, 2, U)
	end
	H0 = copy(H)
	for (m, (ω, α)) in enumerate(bathspecs)
		b = nbands + m
		imp = ((m - 1) % nbands) + 1   # bath mode m couples to this impurity band
		ed_density!(H, b, ω)
		ed_density!(H0, b, ω)
		v = sqrt(α)
		ed_tunneling!(H, b, imp, v)   # a† b
		ed_tunneling!(H, imp, b, v)   # b† a
	end
	cs = fermion_annihilators(nmodes)
	a, adag = cs[1], cs[1]'
	return H, a, adag, H0
end

# Kanamori (multi-orbital SK) impurity with one single-mode bath per band,
# following the KanamoriIM docstring convention: orbital `a` occupies bands
# (2a-1, 2a); bath modes are 2norb+1..4norb.
function kanamori_ed(; U::Real, J::Real, norb::Int, μ::Real, ω::Real=1.0, α::Real=0.5)
	nimp = 2 * norb
	nmodes = 2 * nimp
	dim = 1 << nmodes
	H = zeros(Float64, dim, dim)
	# chemical potential
	for band in 1:nimp
		ed_density!(H, band, μ)
	end
	# intra-orbital interaction
	for a in 1:norb
		i, j = 2*a-1, 2*a
		ed_density_interaction!(H, i, j, U)
	end
	# inter-orbital density-density interaction
	for a in 1:norb, b in 1:norb
		(a != b) && ed_density_interaction!(H, 2*a-1, 2*b, U - 2*J)
	end
	for b in 1:norb, a in (b+1):norb
		ed_density_interaction!(H, 2*a-1, 2*b-1, U - 3*J)
		ed_density_interaction!(H, 2*a,   2*b,   U - 3*J)
	end
	# spin-flip and pair-hopping terms
	for a in 1:norb, b in 1:norb
		if a != b
			ed_interaction!(H, 2*a-1, 2*a, 2*b-1, 2*b, -J)
			ed_interaction!(H, 2*a-1, 2*b, 2*b-1, 2*a, -J)
		end
	end
	H0 = copy(H)
	# bath modes and hybridization
	for m in 1:nimp
		b = nimp + m
		ed_density!(H, b, ω)
		ed_density!(H0, b, ω)
		v = sqrt(α)
		ed_tunneling!(H, b, m, v)
		ed_tunneling!(H, m, b, v)
	end
	cs = fermion_annihilators(nmodes)
	a, adag = cs[1], cs[1]'
	return H, a, adag, H0
end

# BCS impurity model (2 impurity orbitals + single-mode BCS bath), used for
# the bcs code path:
# H = -ϵ_d(n̂₁+n̂₂) + Un̂₁n̂₂ + √α(â₁†ĉ₁ + ĉ₁†â₁ + â₂†ĉ₂ + ĉ₂†â₂)
#     + ω₀(m̂₁+m̂₂) - (Δĉ₁†ĉ₂† + Δ'ĉ₂ĉ₁)
function bcs_ed(U, ϵ_d; ω₀=1, α=0.5, Δ=0.3)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
	Is = one(n̂)
	n_ud = kron(n̂, Is) + kron(Is, n̂)
	nn = kron(n̂,n̂)
	I_ud = kron(Is, Is)
	Himpbare = -ϵ_d*n_ud + U * nn
	Himp = kron(Himpbare, I_ud)
	Hbathbare = ω₀ * n_ud - (Δ * kron(JW*σ₊, σ₊)- conj(Δ) * kron(JW*σ₋, σ₋))
	Hbath = kron(I_ud, Hbathbare)
	tmp = sqrt(α) * (kron(kron(kron(JW*σ₊, JW), σ₋), Is) + kron(Is, kron(JW*σ₊, kron(JW, σ₋))))
	Hhyb = tmp + tmp'
	H = Himp + Hhyb + Hbath
	A, B = kron(kron(σ₋, Is), I_ud), kron(kron(σ₊, Is), I_ud)
	return H, A, B, Himp + Hbath
end

# ----------------------------------------- independent bosons ED builders
# single-particle density matrix n_i for mode i within nbands modes (bit basis)
function ed_density_matrix(nbands::Int, i::Int)
	dim = 1 << nbands
	n = zeros(Float64, dim, dim)
	for s in 0:dim-1
		((s >> (i-1)) & 1) == 1 && (n[s+1, s+1] = 1)
	end
	return n
end

# Fermionic impurity (nbands bands) coupled to a single phonon mode:
#   H = μ Σ n_i + U n₁n₂ + ω₀ b†b + √α (Σ n_i) (b + b†)
# nbands = 1 for U = 0, 2 for U ≠ 0.
function phonon_ed(; μ::Real, U::Real, ω₀::Real=1.0, α::Real=0.5, d::Int=8)
	nbands = (U == zero(U)) ? 1 : 2
	dimf = 1 << nbands
	# fermionic part (impurity + interactions) in the bit basis
	Hf = zeros(Float64, dimf, dimf)
	H0f = zeros(Float64, dimf, dimf)
	for i in 1:nbands
		ed_density!(Hf, i, μ)
		ed_density!(H0f, i, μ)
	end
	(nbands == 2) && (U != zero(U)) && begin
		ed_density_interaction!(Hf, 1, 2, U)
		ed_density_interaction!(H0f, 1, 2, U)
	end
	# bosonic mode
	_b = boson(d=d)
	b̂, b̂′, n̂b = _b["a"], _b["adag"], _b["n"]
	Ib = one(b̂)
	If = one(Hf)
	ntot = zeros(Float64, dimf, dimf)
	for i in 1:nbands
		ntot .+= ed_density_matrix(nbands, i)
	end
	H = kron(Hf, Ib) + ω₀ * kron(If, n̂b) + sqrt(α) * kron(ntot, b̂ + b̂′)
	H0 = kron(H0f, Ib) + ω₀ * kron(If, n̂b)
	cs = fermion_annihilators(nbands)
	A, B = kron(cs[1], Ib), kron(cs[1]', Ib)
	return H, A, B, H0
end

# Fermionic impurity (nbands bands) coupled to BOTH a fermionic single-mode
# bath (one mode per band, (ω_f, α_f)) AND a single phonon mode (ω₀, α):
#   H = μ Σ n_i + U n₁n₂ + Σ_m [ω b†b + √α_f (a†b + h.c.)]
#       + ω₀ c†c + √α (Σ n_i)(c + c†)
function mixedbath_ed(; μ::Real, U::Real, ωf::Real=1.0, αf::Real=0.5, ω₀::Real=1.0, α::Real=0.5, d::Int=8)
	nbands = (U == zero(U)) ? 1 : 2
	nfmodes = 2 * nbands          # impurity modes + fermionic bath modes
	dimf = 1 << nfmodes
	Hf = zeros(Float64, dimf, dimf)
	H0f = zeros(Float64, dimf, dimf)
	for i in 1:nbands
		ed_density!(Hf, i, μ)
		ed_density!(H0f, i, μ)
	end
	(nbands == 2) && (U != zero(U)) && begin
		ed_density_interaction!(Hf, 1, 2, U)
		ed_density_interaction!(H0f, 1, 2, U)
	end
	for m in 1:nbands
		b = nbands + m
		ed_density!(Hf, b, ωf)
		ed_density!(H0f, b, ωf)
		v = sqrt(αf)
		ed_tunneling!(Hf, b, m, v)
		ed_tunneling!(Hf, m, b, v)
	end
	_b = boson(d=d)
	b̂, b̂′, n̂b = _b["a"], _b["adag"], _b["n"]
	Ib = one(b̂)
	If = one(Hf)
	ntot = zeros(Float64, 1 << nfmodes, 1 << nfmodes)
	for i in 1:nbands
		ntot .+= ed_density_matrix(nfmodes, i)
	end
	H = kron(Hf, Ib) + ω₀ * kron(If, n̂b) + sqrt(α) * kron(ntot, b̂ + b̂′)
	H0 = kron(H0f, Ib) + ω₀ * kron(If, n̂b)
	cs = fermion_annihilators(nfmodes)
	A, B = kron(cs[1], Ib), kron(cs[1]', Ib)
	return H, A, B, H0
end
