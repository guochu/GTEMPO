# TDVP-based influence functional construction (fermionic/Grassmann version).
#
# TDVPIF is an influence functional construction algorithm on the same footing
# as XTRGIF and PartialIF. It views the influence functional
# as the "equilibrium state" IF = exp(H) of the influence operator H (the MPO
# form returned by `influenceoperators`), and computes it with a second-order
# single-site TDVP imaginary-time flow
#
#     dz/dτ = H·z ,   τ : 0 → 1 ,
#
# starting from the identity influence functional z(0) = I (β = 0). Each flow
# step is one forward-backward TDVP sweep: in a left-to-right (right-to-left)
# sweep the center tensor AC is evolved by +δτ/2 through Krylov exponentiation
# of the local effective map, factorized by QR (LQ), and the bond matrix C is
# evolved by -δτ/2; the last (first) site performs a full +δτ step.
#
# The initial identity state is zero-padded up to the bond dimension D and
# canonicalized without truncation, so that the zero-weight directions become
# orthonormal directions of the environments and the sweeps can populate the
# full bond profile min(d^j, d^{L-j}, D) as the flow builds up correlations.
# Truncation is applied only in the final canonicalization.
#
# All fermionic signs of the tangent-space effective maps are inherited from
# the `mult` machinery: the AC map reuses `get_left_xy` and the environment
# updates reuse `updatemultleft`/`updatemultright` (with the operator H in the
# "x" slot and the state z in both the bra and "y" slots), and the bond-matrix
# map is realized as the projection AL†·(H·z) of the AC map, so that only
# sign-safe bosonic contractions (bond absorptions and bra-ket overlaps
# analogous to those in `iterativemult`) are added on top.


# the influence operator H driving the flow, in fused-leg GrassmannMPS form:
# `influenceoperators` returns the MPO; multiplying it onto the identity
# influence functional (vacuumstate) fuses its legs into the MPO-as-MPS form
# used by the flow. On lattices whose ordering is not time-local the operator
# is built on the canonical ordering and rotated back with `changeordering`,
# exactly as in `influenceoperatorstepper`.
function _tdvpif_hamiltonian(lattice::ImagGrassmannLattice1Order{O}, corr::ImagCorrelationFunction, alg::TDVPIF; band::Int=1) where O
	T = scalartype(corr)
	if LayoutStyle(lattice) isa TimeLocalLayout
		H = only(influenceoperators(lattice, corr, band=band, algexpan=alg.algexpan))
		return H * vacuumstate(T, lattice)
	else
		lattice2 = similar(lattice, ordering=A1Ā1B1B̄1())
		H = only(influenceoperators(lattice2, corr, band=band, algexpan=alg.algexpan)) * vacuumstate(T, lattice2)
		_, H2 = changeordering(O, lattice2, H, trunc=alg.trunc)
		return H2
	end
end

# on real-time lattices `influenceoperators` returns 4 branch MPOs
# ((+,+), (+,−), (−,+), (−,−)); the total influence operator driving the
# flow is their SUM (MPS/MPO direct sum, bond dimensions add up). Indeed the
# differential IF built by `influenceoperatorstepper` is the Hadamard
# (site-wise) product e^{dt·h1}∘e^{dt·h2}∘e^{dt·h3}∘e^{dt·h4}, and in the
# element-wise algebra of ADT/PT products e^a∘e^b = e^{a+b}, so the generator
# of the full IF is h1+h2+h3+h4 (NOT the Hadamard product h1∘h2∘h3∘h4).
# NOTE: the direct-sum bond dimensions of the four branches add up, so the
# branches are summed one by one, compressing with SVD canonicalization after
# each addition to keep the bond dimension bounded. The compression uses the
# tight `DefaultMPOTruncation` (not `alg.trunc`): the compression error of H
# is exponentially amplified by the flow (IF = e^H), so a loose tolerance
# would degrade the accuracy of the influence functional.
function _tdvpif_hamiltonian(lattice::RealGrassmannLattice{O}, corr::RealCorrelationFunction, alg::TDVPIF; band::Int=1) where O
	T = scalartype(corr)
	if !(OrderingStyle(lattice) isa _AllowedRealGrassmannOrdering)
		lattice2 = similar(lattice, ordering=A1Ā1a1ā1B1B̄1b1b̄1())
		H = _tdvpif_hamiltonian_timelocal(lattice2, corr, alg, T, band=band)
		_, H2 = changeordering(O, lattice2, H, trunc=alg.trunc)
		return H2
	else
		return _tdvpif_hamiltonian_timelocal(lattice, corr, alg, T, band=band)
	end
end

function _tdvpif_hamiltonian_timelocal(lattice::RealGrassmannLattice{<:_AllowedRealGrassmannOrdering}, corr::RealCorrelationFunction, alg::TDVPIF, T; band::Int=1)
	h1, h2, h3, h4 = influenceoperators(lattice, corr, band=band, algexpan=alg.algexpan)
	orth = Orthogonalize(SVD(), DefaultMPOTruncation; normalize=false)
	H = h1 * vacuumstate(T, lattice)
	H = H + h2 * vacuumstate(T, lattice)
	canonicalize!(H, alg=orth)
	H = H + h3 * vacuumstate(T, lattice)
	canonicalize!(H, alg=orth)
	H = H + h4 * vacuumstate(T, lattice)
	canonicalize!(H, alg=orth)
	return H
end

# run the TDVP flow z(τ=1) = e^H·z(0) directly on the input state z, which may
# be the identity influence functional (β = 0) or, more generally, any
# GrassmannMPS on the same lattice, e.g. the impurity dynamics obtained from
# `sysdynamics`; the influence operator H is thereby merged into the impurity
# dynamics in a single flow.
#
# preparation: lift z to the flow bond dimension (zero-padded) and canonicalize
# without truncation, so that the zero-weight directions become orthonormal
# directions of the environments and the sweeps can populate the full bond
# profile min(d^j, d^{L-j}, D); finalization: canonicalize with the truncation
# scheme. The global scaling factor of z is carried through the flow by the
# `_renormalize!` bookkeeping, so the output value is e^H·z(0) regardless of
# the gauge of the input.
function _tdvpif_hybriddynamics!(z::GrassmannMPS, H::GrassmannMPS, alg::TDVPIF)
	increase_bond!(z, alg.trunc.D)
	canonicalize!(z, alg=Orthogonalize(SVD(), NoTruncation(); normalize=false))
	_tdvpif_flow!(z, H, alg)
	# after the flow, a final canonicalization sweep compresses z with SVD
	# truncation to the target bond dimension
	canonicalize!(z, alg=Orthogonalize(SVD(), alg.trunc; normalize=false))
	alg.callback(Float64[])
	return z
end

# absorb the global scaling factor of a GrassmannMPS into its site tensors
# (value = scaling^L · ∏tensors  →  value = 1^L · ∏(scaling·tensors)) and reset
# the scaling to 1, so that downstream code contracting the raw site tensors
# represents the same operator regardless of the gauge
function _absorb_scaling!(x::GrassmannMPS)
	sca = scaling(x)
	(sca == 1) && return x
	for i in 1:length(x)
		x[i] = sca * x[i]
	end
	setscaling!(x, 1)
	return x
end

# effective map on the center tensor at site j: the tangent-space projection
# of the product H·z, built from the left environment hleft::(bra, H, ket),
# the influence operator H[j] and the right environment hright::(ket, H, bra)
function _tdvpif_ac_prime(AC::MPSTensor, Hj::MPSTensor, hleft::MPSTensor, hright::MPSTensor)
	left_xy = get_left_xy(hleft, Hj, AC)
	@tensor mpsj[1,2;5] := left_xy[1,2,3,4] * hright[4,3,5]
	return mpsj
end

# effective map on the bond matrix C in the left-to-right sweep: the projection
# AL†·(H·z) of the AC map at the factorized site, obtained by reforming
# AC = AL·C, applying the AC map and contracting the bra AL back in. All
# contractions are bond absorptions / bra-ket overlaps, which are sign-free.
function _tdvpif_c_prime_left(C::MPSBondTensor, AL::MPSTensor, Hj::MPSTensor, hleft::MPSTensor, hright::MPSTensor)
	AC = @tensor tmp[1,2;4] := AL[1,2,3] * C[3,4]
	ac′ = _tdvpif_ac_prime(AC, Hj, hleft, hright)
	return @tensor c′[3;4] := conj(AL[1,2,3]) * ac′[1,2,4]
end

# effective map on the bond matrix C in the right-to-left sweep: the
# projection AR†·(H·z) of the AC map, with AC = C·AR (C on the left bond)
function _tdvpif_c_prime_right(C::MPSBondTensor, AR::MPSTensor, Hj::MPSTensor, hleft::MPSTensor, hright::MPSTensor)
	AC = @tensor tmp[1,2;4] := C[1,3] * AR[3,2,4]
	ac′ = _tdvpif_ac_prime(AC, Hj, hleft, hright)
	return @tensor c′[1;4] := conj(AR[4,2,3]) * ac′[1,2,3]
end

# initialize the right environments ⟨z|H|z⟩; hstorage[i]::(ket, H, bra) is the
# partial contraction over sites i:L, exactly as in `mult_cache`
function _tdvpif_init_hstorage(z::GrassmannMPS, H::GrassmannMPS)
	L = length(z)
	T = promote_type(scalartype(z), scalartype(H))
	right = ones(T, space_r(z)' ⊗ space_r(H)', space_r(z)')
	hstorage = Vector{typeof(right)}(undef, L+1)
	hstorage[L+1] = right
	hip1 = hstorage[L+1]
	for i in L:-1:2
		hip1 = updatemultright(hip1, z[i], H[i], z[i])
		hstorage[i] = hip1
	end
	hstorage[1] = ones(T, space_l(z) ⊗ space_l(H)', space_l(z))
	return hstorage
end

function _tdvpif_leftsweep!(z::GrassmannMPS, H::GrassmannMPS, hstorage, δτ::Float64, alg::TDVPIF)
	L = length(z)
	krylov = Arnoldi()
	for site in 1:L-1
		(alg.verbosity > 3) && println("TDVPIF: left sweep at site $site")
		# forward half-step of the center tensor
		AC, info = exponentiate(x -> _tdvpif_ac_prime(x, H[site], hstorage[site], hstorage[site+1]), δτ/2, z[site], krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $site"
		AL, C = leftorth!(AC, alg = QR())
		z[site] = AL
		_renormalize!(z, C, false)
		# left environment with the new AL on both the bra and ket sides
		hnew = updatemultleft(hstorage[site], AL, H[site], AL)
		# backward half-step of the bond matrix
		C, info = exponentiate(x -> _tdvpif_c_prime_left(x, AL, H[site], hstorage[site], hstorage[site+1]), -δτ/2, C, krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at bond $(site+1)"
		_renormalize!(z, C, false)
		hstorage[site+1] = hnew
		# absorb the bond matrix into the next site
		z[site+1] = @tensor tmp[1,2;3] := C[1,4] * z[site+1][4,2,3]
	end
	# full step at the last site
	AC, info = exponentiate(x -> _tdvpif_ac_prime(x, H[L], hstorage[L], hstorage[L+1]), δτ, z[L], krylov)
	(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $L"
	z[L] = AC
	_renormalize!(z, z[L], false)
	return z
end

function _tdvpif_rightsweep!(z::GrassmannMPS, H::GrassmannMPS, hstorage, δτ::Float64, alg::TDVPIF)
	krylov = Arnoldi()
	for site in length(z)-1:-1:1
		(alg.verbosity > 3) && println("TDVPIF: right sweep at site $site")
		C, AR = rightorth(z[site+1], (1,), (2, 3), alg=LQ())
		z[site+1] = permute(AR, (1, 2), (3,))
		# right environment with the new AR on both the bra and ket sides
		hnew = updatemultright(hstorage[site+2], z[site+1], H[site+1], z[site+1])
		# backward half-step of the bond matrix; note the projection uses the
		# right environment at bond site+2 (excluding site+1 itself), since the
		# ⟨AR| projection supplies site+1 on both the bra and ket sides
		C, info = exponentiate(x -> _tdvpif_c_prime_right(x, z[site+1], H[site+1], hstorage[site+1], hstorage[site+2]), -δτ/2, C, krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at bond $(site+1)"
		_renormalize!(z, C, false)
		hstorage[site+1] = hnew
		# absorb the bond matrix into the site on the left
		z[site] = @tensor tmp[1,2;4] := z[site][1,2,3] * C[3,4]
		# forward half-step of the center tensor
		AC, info = exponentiate(x -> _tdvpif_ac_prime(x, H[site], hstorage[site], hstorage[site+1]), δτ/2, z[site], krylov)
		(info.converged > 0) || @warn "TDVPIF: Krylov exponentiation failed to converge at site $site"
		z[site] = AC
		_renormalize!(z, z[site], false)
	end
	return z
end

function _tdvpif_flow!(z::GrassmannMPS, H::GrassmannMPS, alg::TDVPIF)
	# The flow contracts the raw site tensors of `H` and never applies its
	# global scaling factor. The GrassmannMPS convention is
	#     value = scaling^L · ∏(site tensors),
	# so any H with scaling ≠ 1 (e.g. after a canonicalization, which
	# redistributes local weights into the scaling factor) would silently be
	# evolved as H / scaling^L. Absorb the scaling into the site tensors first
	# to make the flow independent of the gauge in which H is represented.
	_absorb_scaling!(H)
	hstorage = _tdvpif_init_hstorage(z, H)
	nsteps = round(Int, 1 / alg.δ)
	δτ = 1 / nsteps
	for n in 1:nsteps
		_tdvpif_leftsweep!(z, H, hstorage, δτ, alg)
		_tdvpif_rightsweep!(z, H, hstorage, δτ, alg)
		(alg.verbosity > 1) && println("TDVPIF step $n/$nsteps, τ = $(n * δτ)")
	end
	(alg.verbosity > 1) && println("TDVPIF flow finished, τ = 1")
	return z
end

"""
	hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::TDVPIF; band::Int=1)

Construct the influence functional with the `TDVPIF` algorithm: evolve the identity influence functional along the imaginary-time flow dz/dτ = H·z (H = influence operator) from τ = 0 to τ = 1 with second-order single-site TDVP sweeps.

# Returns
The influence functional, represented as a `GrassmannMPS`.

See also [`TDVPIF`](@ref).
"""
function hybriddynamics(lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::TDVPIF; band::Int=1)
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	T = promote_type(scalartype(corr), scalartype(lattice), Float64)
	z = vacuumstate(T, lattice)
	return hybriddynamics!(z, lattice, corr, alg, band=band)
end

"""
	hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::TDVPIF; band::Int=1)

In-place version of the `TDVPIF` algorithm: the influence operator H drives the TDVP flow directly on `gmps`, i.e. the flow evolves `z(τ=1) = e^H·gmps` from `z(0) = gmps` in place. This allows merging the influence functional into an arbitrary initial GrassmannMPS — e.g. the impurity dynamics obtained from `sysdynamics` — in a single flow instead of constructing the IF separately and multiplying it in afterwards.

# Returns
The modified `gmps`.

See also [`TDVPIF`](@ref).
"""
function hybriddynamics!(gmps::GrassmannMPS, lattice::AbstractGrassmannLattice, corr::AbstractCorrelationFunction, alg::TDVPIF; band::Int=1)
	(1 <= band <= lattice.bands) || throw(BoundsError(1:lattice.bands, band))
	H = _tdvpif_hamiltonian(lattice, corr, alg, band=band)
	return _tdvpif_hybriddynamics!(gmps, H, alg)
end
