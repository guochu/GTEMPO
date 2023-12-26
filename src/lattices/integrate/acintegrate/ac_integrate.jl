
function _ac_integrate(alg::ExactIntegrate, lattice::AbstractGrassmannLattice, x::GrassmannMPS...; trunc::TruncationScheme=DefaultIntegrationTruncation)
	(ConjugationStyle(lattice) isa AdjacentConjugation) || throw(ArgumentError("AdjacentConjugation expected"))
	(all(v->length(v)==length(lattice), x)) || throw(DimensionMismatch())
	Lhalf = div(length(lattice), 2)
	tmp2 = _l_LL(x...)
	for i in 1:Lhalf
		tmp2 = update_pair_left(tmp2, i, x..., trunc=trunc)
	end
	return TK.scalar(tmp2)
end


pos2pairindex(pos::Int) = div(pos+1, 2)
