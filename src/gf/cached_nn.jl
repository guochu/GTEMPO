
function cached_nn(lattice::AbstractGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                    cache::AbstractExpectationCache=environments(lattice, A, B...), b1::Symbol=:τ, b2::Symbol=b1, 
                    band::Int=1, kwargs...)
	a1, a2 = get_nn_contour_pos(lattice, i, band, b1)
	if i == j
		return cached_gf(lattice, (a1, a2), A, B...; cache=cache, kwargs...)
	end
	a3, a4 = get_nn_contour_pos(lattice, j, band, b2)
	return cached_gf(lattice, (a1, a2, a3, a4), A, B...; cache=cache, kwargs...)
end
