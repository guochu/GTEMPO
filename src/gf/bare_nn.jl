
"""
	nn(lattice::AbstractGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
		b1::Symbol=:τ, b2::Symbol=:τ, band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), 
		Z::Real = integrate(lattice, A, B..., alg=alg))

Calculate the density-density correlation
"""
function nn(lattice::AbstractGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
			b1::Symbol=:τ, b2::Symbol=b1, band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
	a1, a2 = get_nn_contour_pos(lattice, i, band, b1)
	if i == j
		return gf(lattice, (a1, a2), A, B...; alg=alg, Z=Z)
	end
	a3, a4 = get_nn_contour_pos(lattice, j, band, b2)
    return gf(lattice, (a1, a2, a3, a4), A, B...; alg=alg, Z=Z)
end



function get_nn_contour_pos(lattice::AbstractGrassmannLattice, j::Int, band::Int, b::Symbol)
	if b == :-
		return ContourIndex(j, conj=true, branch=b, band=band), ContourIndex(j+1, conj=false, branch=b, band=band)
	else
		return ContourIndex(j+1, conj=true, branch=b, band=band), ContourIndex(j, conj=false, branch=b, band=band)
	end
end