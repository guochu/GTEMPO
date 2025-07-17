
function _ac_integrate(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, x::GrassmannMPS...)
	(ConjugationStyle(lattice) isa AdjacentConjugation) || throw(ArgumentError("AdjacentConjugation expected"))
	transfer = GrassmannTransferMatrix(x...)
	tmp2 = l_LL(x...) * transfer
	return TK.scalar(tmp2)
end


pos2pairindex(pos::Int) = div(pos+1, 2)

