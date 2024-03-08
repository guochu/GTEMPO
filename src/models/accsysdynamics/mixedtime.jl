# is there a point to implement this? For SIAM we have already got the exact expression,
# but for more orbitals the computational cost would be too high

function accsysdynamics_fast(lattice::MixedGrassmannLattice{O}, model::AbstractImpurityHamiltonian; scaling::Int=10, trunc::TruncationScheme=DefaultKTruncation) where O
	error("accsysdynamics_fast not implemented for MixedGrassmannLattice")
end
