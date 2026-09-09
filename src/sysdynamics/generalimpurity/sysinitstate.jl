# Initial states of the impurity
# -------------------------------
# The impurity initial state (an operator in the Fock basis, e.g. a
# density matrix) is converted directly into a SparseGMPS on the lattice
# (see `_tosparsegmps`) and multiplied into the target GMPS with the
# sparse mult!. The bra variables sit on the forward branch and the ket
# variables on the backward branch at time slice 1.

"""
	sysinitialstate!(gmps, lattice, fockstate::FockMatrix; trunc) -> GrassmannMPS

In-place: multiply the Fock-space operator `fockstate` (e.g. a density
matrix in the occupation number basis) into the GrassmannMPS `gmps`. The
operator is converted directly into a SparseGMPS (see `_tosparsegmps`)
with the bra (conjugated) variables on the forward branch and the ket
variables on the backward branch at time slice 1, and multiplied into
`gmps` with the sparse `mult!` using the truncation scheme `trunc`.
"""
function sysinitialstate!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, fockstate::FockMatrix;
						trunc::TruncationScheme=DefaultKTruncation)
    (fockstate.bands == lattice.bands) || throw(DimensionMismatch("FockMatrix bands $(fockstate.bands) do not match lattice bands $(lattice.bands)"))
    M = lattice.bands
    bpos = [index(lattice, 1, conj=true, branch=:+, band=i) for i in 1:M]
    kpos = [index(lattice, 1, conj=false, branch=:-, band=i) for i in 1:M]
    sparse = _tosparsegmps(lattice, fockstate, bpos, kpos)
    return mult!(gmps, sparse, trunc=trunc)
end

"""
	sysinitialstate(lattice, fockstate::FockMatrix; trunc) -> GrassmannMPS

The Fock-space operator `fockstate` as a GrassmannMPS on `lattice`,
built by applying `sysinitialstate!` to the vacuum state.
"""
sysinitialstate(lattice::RealGrassmannLattice, fockstate::FockMatrix; kwargs...) =
	sysinitialstate!(vacuumstate(lattice), lattice, fockstate; kwargs...)

"""
	systhermalstate!(gmps, lattice, model; β, trunc) -> GrassmannMPS

In-place: multiply the normalized impurity thermal equilibrium state
`exp(-βĤ)/tr(exp(-βĤ))` (for `β == Inf` the ground state projector, see
`fock_thermalstate`) into the GrassmannMPS `gmps`.
"""
function systhermalstate!(gmps::GrassmannMPS, lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian;
							β::Real, trunc::TruncationScheme=DefaultKTruncation)
    return sysinitialstate!(gmps, lattice, fock_thermalstate(model, β); trunc=trunc)
end

"""
	systhermalstate(lattice, model; β, trunc) -> GrassmannMPS

The impurity thermal equilibrium state as a GrassmannMPS on `lattice`,
built by applying `systhermalstate!` to the vacuum state.
"""
systhermalstate(lattice::RealGrassmannLattice, model::AbstractImpurityHamiltonian; kwargs...) =
	systhermalstate!(vacuumstate(lattice), lattice, model; kwargs...)
