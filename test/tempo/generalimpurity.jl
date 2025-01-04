println("------------------------------------")
println("|          General Impurity        |")
println("------------------------------------")

function gSISB(;μ::Real, U::Real)
	if U != zero(U)
		h = ImpurityHamiltonian(bands=2)
		push!(h, interaction(1,2,2,1, coeff=U))
	else
		h = ImpurityHamiltonian(bands=1)
	end
	for band in 1:h.bands
		push!(h, tunneling(band, band, coeff=μ))
	end
	return h
end

function gSKIM(; U::Real, J::Real, norb::Int, μ::Real=-U/2)
	h = ImpurityHamiltonian(bands = 2*norb)

	for band in 1:h.bands
		push!(h, tunneling(band, band, coeff=μ))
	end
	for a in 1:norb
		i, j = 2 * a - 1, 2 * a
		push!(h, interaction(i,j,j,i, coeff=U))
	end
	for a in 1:norb
		for b in 1:norb
			if a != b
				i, j = 2*a-1, 2*b
				push!(h, interaction(i,j,j,i, coeff=U-2*J))
			end
		end
	end
	for b in 1:norb
		for a in (b+1):norb
			i, j = 2*a-1, 2*b-1
			push!(h, interaction(i,j,j,i, coeff=U-3*J))
			i, j = 2*a, 2*b
			push!(h, interaction(i,j,j,i, coeff=U-3*J))
		end
	end
	for a in 1:norb
		for b in 1:norb
			if a != b
				i, j, k, l = 2*a-1, 2*a, 2*b-1, 2*b
				push!(h, interaction(i,j,k,l, coeff=-J))
				i, j, k, l = 2*a-1, 2*b, 2*b-1, 2*a
				push!(h, interaction(i,j,k,l, coeff=-J))
			end
		end
	end
	return h
end

@testset "SIAM: imag-time" begin
	ϵ_d = 0.7
	δτ = 0.1
	β = 1
	N = round(Int, β / δτ)
	bath = fermionicbath(spectrum_func(1), β=β, μ=0.5)

	rtol = 1.0e-3
	tol = 1.0e-4

	# 1 orb
	for U in [0, 1]
		model1 = SISB(bath, U=U, μ=ϵ_d)
		model2 = gSISB(U=U, μ=ϵ_d)
		bands = (U == zero(U)) ? 1 : 2
		lattice = GrassmannLattice(δτ=δτ, N=N, bands=bands, contour=:imag)
		K1 = accsysdynamics(lattice, model1)
		K2 = accsysdynamics_fast(lattice, model2, scaling=1000)
		@test distance(K1, K2) / norm(K1) < rtol
		for band in 1:lattice.bands
			K1 = boundarycondition!(K1, lattice, band=band)
			K2 = boundarycondition!(K2, lattice, band=band)
		end
		cache1 = environments(lattice, K1)
		cache2 = environments(lattice, K2)
		for band in 1:lattice.bands
			g1 = cached_Gτ(lattice, K1, cache=cache1, band=band)
			g2 = cached_Gτ(lattice, K2, cache=cache2, band=band)
			@test norm(g1-g2) / norm(g1) < tol
		end
	end
end

@testset "SIAM: real-time" begin
	ϵ_d = 0.5
	δt = 0.1
	N = 5
	bath = fermionicbath(spectrum_func(1), β=10, μ=0.5)

	rtol = 1.0e-2
	tol = 1.0e-2

	# 1 orb
	for U in [0, 1]
		model1 = SISB(bath, U=U, μ=ϵ_d)
		model2 = gSISB(U=U, μ=ϵ_d)
		bands = (U == zero(U)) ? 1 : 2
		lattice = GrassmannLattice(δt=δt, N=N, bands=bands, contour=:real)
		K1 = accsysdynamics(lattice, model1)
		K2 = accsysdynamics_fast(lattice, model2, scaling=100)
		@test distance(K1, K2) / norm(K1) < rtol
		for band in 1:lattice.bands
			K1 = boundarycondition!(K1, lattice, band=band)
			K2 = boundarycondition!(K2, lattice, band=band)
		end
		cache1 = environments(lattice, K1)
		cache2 = environments(lattice, K2)
		for band in 1:lattice.bands
			for i in 1:lattice.k, j in 1:lattice.k
				for b1 in (:+, :-), b2 in (:+, :-), c1 in (true, false)
					g1 = cached_Gt(lattice, i, j, K1, b1=b1, b2=b2, c1=c1, c2=!c1, cache=cache1, band=band)
					g2 = cached_Gt(lattice, i, j, K2, b1=b1, b2=b2, c1=c1, c2=!c1, cache=cache2, band=band)
					@test _error(g1, g2, tol) < tol
				end
			end
		end
	end
end

@testset "SKIM: imag-time" begin
	U = 1.
	ϵ_d = 0.7
	J = 1.1

	δτ = 0.1
	β = 1
	N = round(Int, β / δτ)
	bath = fermionicbath(spectrum_func(1), β=β, μ=0.5)

	rtol = 1.0e-2
	tol = 1.0e-2

	for norb in [1, 2]
		model1 = SKIM(bath, U=U, μ=ϵ_d, J=J, norb=norb)
		model2 = gSKIM(U=U, μ=ϵ_d, J=J, norb=norb)
		bands = 2 * norb
		lattice = GrassmannLattice(δτ=δτ, N=N, bands=bands, contour=:imag)
		K1 = accsysdynamics(lattice, model1)
		K2 = accsysdynamics(lattice, model2)
		@test distance(K1, K2) / norm(K1) < rtol
		for band in 1:lattice.bands
			K1 = boundarycondition!(K1, lattice, band=band)
			K2 = boundarycondition!(K2, lattice, band=band)
		end
		cache1 = environments(lattice, K1)
		cache2 = environments(lattice, K2)
		for band in 1:lattice.bands
			g1 = cached_Gτ(lattice, K1, cache=cache1, band=band)
			g2 = cached_Gτ(lattice, K2, cache=cache2, band=band)
			@test norm(g1-g2) / norm(g1) < tol
		end
	end
	
end

@testset "SKIM: real-time" begin
	U = 1.
	ϵ_d = 0.7
	J = 1.1

	δt = 0.1
	N = 3
	bath = fermionicbath(spectrum_func(1), β=10, μ=0.5)

	rtol = 1.0e-2
	tol = 1.0e-2

	# 1 orb
	for norb in [1]
		model1 = SKIM(bath, U=U, μ=ϵ_d, J=J, norb=norb)
		model2 = gSKIM(U=U, μ=ϵ_d, J=J, norb=norb)
		bands = 2 * norb
		lattice = GrassmannLattice(δt=δt, N=N, bands=bands, contour=:real)
		K1 = accsysdynamics(lattice, model1)
		K2 = accsysdynamics(lattice, model2)
		@test distance(K1, K2) / norm(K1) < rtol
		for band in 1:lattice.bands
			K1 = boundarycondition!(K1, lattice, band=band)
			K2 = boundarycondition!(K2, lattice, band=band)
		end
		cache1 = environments(lattice, K1)
		cache2 = environments(lattice, K2)
		for band in 1:lattice.bands
			for i in 1:lattice.k, j in 1:lattice.k
				for b1 in (:+, :-), b2 in (:+, :-), c1 in (true, false)
					g1 = cached_Gt(lattice, i, j, K1, b1=b1, b2=b2, c1=c1, c2=!c1, cache=cache1, band=band)
					g2 = cached_Gt(lattice, i, j, K2, b1=b1, b2=b2, c1=c1, c2=!c1, cache=cache2, band=band)
					@test _error(g1, g2, tol) < tol
				end
			end
		end
	end
end

