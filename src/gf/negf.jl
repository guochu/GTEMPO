# calculate all the green's functions, which are used for non-equilibrium DMFT

function cached_Gm(lattice::MixedGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
					cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, kwargs...)
	T = scalartype(lattice)
	# G₁₁
	G11 = zeros(T, lattice.kt, lattice.kt)
	G22 = zeros(T, lattice.kt, lattice.kt)
	G12 = zeros(T, lattice.kt, lattice.kt)
	G21 = zeros(T, lattice.kt, lattice.kt)
	for i in 1:lattice.kt, j in 1:lattice.kt
		G12[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:+, b2=:-, band=band, kwargs...)
		G21[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
	end
	for i in 1:lattice.kt
		for j in 1:i-1
			G11[i, j] = G21[i, j]
		end
		for j in i:lattice.kt
			G11[i, j] = G12[i, j]
		end
		for j in 1:i
			G22[i, j] = G12[i, j]
		end
		for j in (i+1):lattice.kt
			G22[i, j] = G21[i, j]
		end
		
	end
	G13 = zeros(T, lattice.kt, lattice.Nτ)
	G31 = zeros(T, lattice.Nτ, lattice.kt)

	for i in 1:lattice.kt, j in 1:lattice.Nτ
		G13[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:+, b2=:τ, band=band, kwargs...)
		# G23[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:-, b2=:τ, band=band, kwargs...)
	end
	G23 = G13
	for i in 1:lattice.Nτ, j in 1:lattice.kt
		G31[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:τ, b2=:+, band=band, kwargs...)
	end
	G32 = G31
	
	G33 = zeros(T, lattice.Nτ, lattice.Nτ)
	for i in 1:lattice.Nτ, j in 1:lattice.Nτ
		G33[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:τ, b2=:τ, band=band, kwargs...)
	end

	r = Matrix{Matrix{T}}(undef, 3, 3)
	r[1,1] = G11
	r[1,2] = G12
	r[1,3] = G13
	r[2,1] = G21
	r[2,2] = G22
	r[2,3] = G23
	r[3,1] = G31
	r[3,2] = G32
	r[3,3] = G33

	return r
end

function cached_Gm_full(lattice::MixedGrassmannLattice, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
					cache::AbstractExpectationCache=environments(lattice, A, B...), band::Int=1, kwargs...)
	T = scalartype(lattice)
	# G₁₁
	G11 = zeros(T, lattice.kt, lattice.kt)
	G22 = zeros(T, lattice.kt, lattice.kt)
	G12 = zeros(T, lattice.kt, lattice.kt)
	G21 = zeros(T, lattice.kt, lattice.kt)
	for i in 1:lattice.kt, j in 1:lattice.kt
		G11[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:+, b2=:+, band=band, kwargs...)
		G22[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:-, b2=:-, band=band, kwargs...)
		G12[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:+, b2=:-, band=band, kwargs...)
		G21[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:-, b2=:+, band=band, kwargs...)
	end
	G13 = zeros(T, lattice.kt, lattice.kτ)
	G23 = zeros(T, lattice.kt, lattice.kτ)
	G31 = zeros(T, lattice.kτ, lattice.kt)
	G32 = zeros(T, lattice.kτ, lattice.kt)

	for i in 1:lattice.kt, j in 1:lattice.kτ
		G13[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:+, b2=:τ, band=band, kwargs...)
		G23[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:-, b2=:τ, band=band, kwargs...)
	end
	for i in 1:lattice.kτ, j in 1:lattice.kt
		G31[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:τ, b2=:+, band=band, kwargs...)
		G32[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:τ, b2=:-, band=band, kwargs...)		
	end
	G33 = zeros(T, lattice.kτ, lattice.kτ)
	for i in 1:lattice.kτ, j in 1:lattice.kτ
		G33[i, j] = cached_contour_ordered_Gm(lattice, i, j, A, B...; cache=cache, b1=:τ, b2=:τ, band=band, kwargs...)
	end

	r = Matrix{Matrix{T}}(undef, 3, 3)
	r[1,1] = G11
	r[1,2] = G12
	r[1,3] = G13
	r[2,1] = G21
	r[2,2] = G22
	r[2,3] = G23
	r[3,1] = G31
	r[3,2] = G32
	r[3,3] = G33

	return r
end

