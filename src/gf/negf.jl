# calculate all the green's functions, which are used for non-equilibrium DMFT


function cached_GF(lattice::MixedGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
					cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, kwargs...)
	T = scalartype(lattice)
	# G₁₁
	G11 = zeros(T, lattice.kt, lattice.kt)
	G22 = zeros(T, lattice.kt, lattice.kt)
	G12 = zeros(T, lattice.kt, lattice.kt)
	G21 = zeros(T, lattice.kt, lattice.kt)
	for i in 1:lattice.kt, j in 1:lattice.kt
		G11[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:+, b2=:+, band=band, kwargs...)
		G22[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)
		G12[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:+, b2=:-, band=band, kwargs...)
		G21[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
	end
	G13 = zeros(T, lattice.kt, lattice.kτ)
	G23 = zeros(T, lattice.kt, lattice.kτ)
	G31 = zeros(T, lattice.kτ, lattice.kt)
	G32 = zeros(T, lattice.kτ, lattice.kt)

	for i in 1:lattice.kt, j in 1:lattice.kτ
		G13[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:+, b2=:τ, band=band, kwargs...)
		G23[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:-, b2=:τ, band=band, kwargs...)
		G31[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:τ, b2=:+, band=band, kwargs...)
		G32[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:τ, b2=:-, band=band, kwargs...)
	end
	G33 = zeros(T, lattice.kτ, lattice.kτ)
	for i in 1:lattice.kτ, j in 1:lattice.kτ
		G33[i, j] = cached_contour_ordered_gf(lattice, i, j, A, B...; cache=cache, b1=:τ, b2=:τ, band=band, kwargs...)
	end

	return [G11 G12 G13; G21 G22 G23; G31 G32 G33]
end

