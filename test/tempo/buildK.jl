println("------------------------------------")
println("|      Accurate GMPS for K         |")
println("------------------------------------")

function _error(a, b, tol)
	ab = a - b
	if abs(a) < tol
		return abs(ab)
	else
		return abs(ab / a)
	end
end

@testset "change Grassmann ordering" begin
	U = 1.
	ϵ_d = 0.3
	δτ = 0.02
	β = 0.06
	N = round(Int, β / δτ)
	bath = fermionicbath(spectrum_func(1), β=β, μ=1)
	tol = 1.0e-6

	for norb in 1:2
		exact_model = SKIM(bath, U=U, μ=ϵ_d, J=1.1, norb=norb)
		bands = 2 * norb
		lattice1 = GrassmannLattice(δt=δτ, N=N, bands=bands, contour=:real, ordering=A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2())
		lattice2 = GrassmannLattice(δt=δτ, N=N, bands=bands, contour=:real, ordering=A1A1a1a1B1B1b1b1())
		for f in (true, false)
			K1 = sysdynamics(lattice1, exact_model, forward=f)
			K2 = sysdynamics(lattice2, exact_model, forward=f)
			K2′ = changeordering(A1A1a1a1B1B1b1b1, lattice1, K1)[2]
			@test distance(K2, K2′) / norm(K2) < tol
		end
		if norb == 1
			K1 = sysdynamics(lattice1, exact_model)
			K2 = sysdynamics(lattice2, exact_model)
			K2′ = changeordering(A1A1a1a1B1B1b1b1, lattice1, K1)[2]
			@test distance(K2, K2′) / norm(K2) < tol			
		end
	end
	

end

@testset "Zoomout grassmannmps" begin
	δτ = 0.02
	N = 5
	tol = 1.0e-10
	for ordering in (A1A1B1B1(), A1B1B1A1())
		for scaling in (2, 3)
			for bands in (1,2,3)
				lattice = GrassmannLattice(δτ=δτ, N=N, bands=bands, contour=:imag, ordering=ordering)
				lattice_scaling = zoomin(lattice, scaling=scaling)
				K1 = randomgmps(scalartype(lattice_scaling), length(lattice_scaling), D=6)
				canonicalize!(K1, alg=Orthogonalize(trunc=NoTruncation()))
				K2 = zoomout(K1, lattice_scaling, scaling=scaling)

				Z1 = integrate(lattice_scaling, K1)
				Z2 = integrate(lattice, K2)
				@test _error(Z1, Z2, tol) < tol
				for i in 1:lattice.N
					for band in 1:lattice.bands
						g1 = gf(lattice_scaling, (i-1)*scaling+1, K1, band=band, Z=Z1)
						g2 = gf(lattice, i, K2, band=band, Z=Z2)
						@test _error(g1, g2, tol) < tol
					end
				end	
			end
		end
	end

	for ordering in (A1A1a1a1B1B1b1b1(), A1a1B1b1b1B1a1A1(), A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2(), A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2())
		for scaling in (2, 3)
			for bands in (1,2,3)
				lattice = GrassmannLattice(δt=δτ, N=N, bands=bands, contour=:real, ordering=ordering)
				lattice_scaling = zoomin(lattice, scaling=scaling)
				K1 = randomgmps(scalartype(lattice_scaling), length(lattice_scaling), D=4)
				canonicalize!(K1, alg=Orthogonalize(trunc=NoTruncation()))
				K2 = zoomout(K1, lattice_scaling, scaling=scaling)

				Z1 = integrate(lattice_scaling, K1)
				Z2 = integrate(lattice, K2)
				@test _error(Z1, Z2, tol) < tol
				for i in 1:lattice.k-1
					for j in 1:lattice.k-1
						for f1 in (true, false), f2 in (true, false)
							for band in 1:lattice.bands
								g1 = gf(lattice_scaling, (i-1)*scaling+1, (j-1)*scaling+1, K1, f1=f1, f2=f2, band=band, Z=Z1)
								g2 = gf(lattice, i, j, K2, f1=f1, f2=f2, band=band, Z=Z2)
								@test _error(g1, g2, tol) < tol
							end		
						end
					end
				end
			end
		end
	end
end

@testset "build K-imaginary time" begin
	U = 1.
	ϵ_d = 0.7
	δτ = 0.02
	β = 0.1
	N = round(Int, β / δτ)
	bath = fermionicbath(spectrum_func(1), β=β, μ=1)
	tol = 1.0e-9

	# 1 orb
	exact_model = SISB(bath, U=U, μ=ϵ_d)

	for ordering in (A1A1B1B1(), A1B1B1A1())
		lattice = GrassmannLattice(δτ=δτ, N=N, bands=2, contour=:imag, ordering=ordering)
		for scaling in [3,4]
			lattice_scaling = zoomin(lattice, scaling=scaling)
			K1 = sysdynamics(lattice_scaling, exact_model)
			for band in 1:lattice.bands
				K1 = boundarycondition(K1, lattice_scaling, band=band)
			end
			K2 = accsysdynamics(lattice, exact_model, scaling=scaling)
			for band in 1:lattice.bands
				K2 = boundarycondition(K2, lattice, band=band)
			end	
			K3 = accsysdynamics_fast(lattice, exact_model, scaling=scaling)
			for band in 1:lattice.bands
				K3 = boundarycondition(K3, lattice, band=band)
			end	
			Z1 = integrate(lattice_scaling, K1)
			Z2 = integrate(lattice, K2)
			@test _error(Z1, Z2, tol) < tol
			Z3 = integrate(lattice, K3)
			@test _error(Z1, Z3, tol) < tol
			for i in 1:lattice.N
				for band in 1:lattice.bands
					g1 = gf(lattice_scaling, (i-1)*scaling+1, K1, band=band, Z=Z1)
					g2 = gf(lattice, i, K2, band=band, Z=Z2)
					@test _error(g1, g2, tol) < tol
					g3 = gf(lattice, i, K3, band=band, Z=Z3)
					@test _error(g1, g3, tol) < tol
				end
			end
		end
	end

	# multi orb
	for norb in 1:2
		exact_model = SKIM(bath, U=U, μ=ϵ_d, J=1.1, norb=norb)
		for ordering in (A1A1B1B1(), A1B1B1A1())
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=2*norb, contour=:imag, ordering=ordering)
			for scaling in [2,3]
				lattice_scaling = zoomin(lattice, scaling=scaling)
				K1 = sysdynamics(lattice_scaling, exact_model)
				for band in 1:lattice.bands
					K1 = boundarycondition(K1, lattice_scaling, band=band)
				end
				K2 = accsysdynamics(lattice, exact_model, scaling=scaling)
				for band in 1:lattice.bands
					K2 = boundarycondition(K2, lattice, band=band)
				end	
				K3 = accsysdynamics_fast(lattice, exact_model, scaling=scaling)
				for band in 1:lattice.bands
					K3 = boundarycondition(K3, lattice, band=band)
				end	
				Z1 = integrate(lattice_scaling, K1)
				Z2 = integrate(lattice, K2)
				@test _error(Z1, Z2, tol) < tol
				Z3 = integrate(lattice, K3)
				@test _error(Z1, Z3, tol) < tol
				for i in 1:lattice.N
					for band in 1:lattice.bands
						g1 = gf(lattice_scaling, (i-1)*scaling+1, K1, band=band, Z=Z1)
						g2 = gf(lattice, i, K2, band=band, Z=Z2)
						@test _error(g1, g2, tol) < tol
						g3 = gf(lattice, i, K3, band=band, Z=Z3)
						@test _error(g1, g3, tol) < tol
					end
				end
			end			
		end
	end
	
end

@testset "build K-real time" begin
	U = 1.
	ϵ_d = 0.3
	δτ = 0.02
	β = 0.06
	N = round(Int, β / δτ)
	bath = fermionicbath(spectrum_func(1), β=β, μ=1)
	tol = 1.0e-9

	# 1 orb, time-local ordering
	exact_model = SISB(bath, U=U, μ=ϵ_d)
	for ordering in (A1A1a1a1B1B1b1b1(), A1a1B1b1b1B1a1A1())
		lattice = GrassmannLattice(δt=δτ, N=N, bands=2, contour=:real, ordering=ordering)
		for scaling in [3,4]
			lattice_scaling = zoomin(lattice, scaling=scaling)
			K1 = sysdynamics(lattice_scaling, exact_model)
			for band in 1:lattice.bands
				K1 = boundarycondition(K1, lattice_scaling, band=band)
			end
			K2 = accsysdynamics(lattice, exact_model, scaling=scaling)
			for band in 1:lattice.bands
				K2 = boundarycondition(K2, lattice, band=band)
			end	
			Z1 = integrate(lattice_scaling, K1)
			Z2 = integrate(lattice, K2)
			@test _error(Z1, Z2, tol) < tol
			for i in 1:lattice.k-1
				for j in 1:lattice.k-1
					for f1 in (true, false), f2 in (true, false)
						for band in 1:lattice.bands
							g1 = gf(lattice_scaling, (i-1)*scaling+1, (j-1)*scaling+1, K1, f1=f1, f2=f2, band=band,Z=Z1)
							g2 = gf(lattice, i, j, K2, f1=f1, f2=f2, band=band,Z=Z2)
							@test _error(g1, g2, tol) < tol
						end		
					end
				end
			end
		end
	end

	# 1 orb, branch-local ordering
	exact_model = SISB(bath, U=U, μ=ϵ_d)
	for ordering in (A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2(), A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2())
		lattice = GrassmannLattice(δt=δτ, N=N, bands=2, contour=:real, ordering=ordering)
		for scaling in [3,4]
			lattice_scaling = zoomin(lattice, scaling=scaling)
			K1 = sysdynamics(lattice_scaling, exact_model)
			for band in 1:lattice.bands
				K1 = boundarycondition(K1, lattice_scaling, band=band)
			end
			K2 = accsysdynamics(lattice, exact_model, scaling=scaling)
			for band in 1:lattice.bands
				K2 = boundarycondition(K2, lattice, band=band)
			end	
			K3 = accsysdynamics_fast(lattice, exact_model, scaling=scaling)
			K3 = [K3]
			for band in 1:lattice.bands
				K3 = boundarycondition_branching(K3, lattice, band=band)
			end	
			Z1 = integrate(lattice_scaling, K1)
			Z2 = integrate(lattice, K2)
			@test _error(Z1, Z2, tol) < tol
			Z3 = integrate(lattice, K3)
			@test _error(Z1, Z3, tol) < tol
			for i in 1:lattice.k-1
				for j in 1:lattice.k-1
					for f1 in (true, false), f2 in (true, false)
						for band in 1:lattice.bands
							g1 = gf(lattice_scaling, (i-1)*scaling+1, (j-1)*scaling+1, K1, f1=f1, f2=f2, band=band,Z=Z1)
							g2 = gf(lattice, i, j, K2, f1=f1, f2=f2, band=band,Z=Z2)
							@test _error(g1, g2, tol) < tol
							g3 = gf(lattice, i, j, K3, f1=f1, f2=f2, band=band,Z=Z3)
							@test _error(g1, g3, tol) < tol
						end		
					end
				end
			end
		end
	end

	for norb in 1:2
		exact_model = SKIM(bath, U=U, μ=ϵ_d, J=1.1, norb=norb)
		for ordering in (A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2(), A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2())
			lattice = GrassmannLattice(δt=δτ, N=N, bands=2*norb, contour=:real, ordering=ordering)
			for scaling in [2]
				lattice_scaling = zoomin(lattice, scaling=scaling)
				K1 = [sysdynamics(lattice_scaling, exact_model)]
				for band in 1:lattice.bands
					K1 = boundarycondition_branching(K1, lattice_scaling, band=band)
				end
				K2 = [accsysdynamics(lattice, exact_model, scaling=scaling)]
				for band in 1:lattice.bands
					K2 = boundarycondition_branching(K2, lattice, band=band)
				end	
				K3 = [accsysdynamics_fast(lattice, exact_model, scaling=scaling)]
				for band in 1:lattice.bands
					K3 = boundarycondition_branching(K3, lattice, band=band)
				end	
				Z1 = integrate(lattice_scaling, K1)
				Z2 = integrate(lattice, K2)
				@test _error(Z1, Z2, tol) < tol
				Z3 = integrate(lattice, K3)
				@test _error(Z1, Z3, tol) < tol
				for i in 1:lattice.k-1
					for j in 1:lattice.k-1
						for f1 in (true, false), f2 in (true, false)
							for band in 1:lattice.bands
								g1 = gf(lattice_scaling, (i-1)*scaling+1, (j-1)*scaling+1, K1, f1=f1, f2=f2, band=band,Z=Z1)
								g2 = gf(lattice, i, j, K2, f1=f1, f2=f2, band=band,Z=Z2)
								@test _error(g1, g2, tol) < tol
								g3 = gf(lattice, i, j, K3, f1=f1, f2=f2, band=band,Z=Z3)
								@test _error(g1, g3, tol) < tol
							end		
						end
					end
				end
			end
		end
	end

	ordering = A1A1a1a1B1B1b1b1()
	for norb in 1:2
		exact_model = SKIM(bath, U=U, μ=ϵ_d, J=1.1, norb=norb)
		lattice = GrassmannLattice(δt=δτ, N=N, bands=2*norb, contour=:real, ordering=ordering)
		for scaling in [2]
			for f in (true, false)
				lattice_scaling = zoomin(lattice, scaling=scaling)
				K1 = sysdynamics(lattice_scaling, exact_model, forward=f)
				for band in 1:lattice.bands
					K1 = boundarycondition(K1, lattice_scaling, band=band)
				end
				K2 = accsysdynamics_fast(lattice, exact_model, scaling=scaling, forward=f)
				for band in 1:lattice.bands
					K2 = boundarycondition(K2, lattice, band=band)
				end	
				Z1 = integrate(lattice_scaling, K1)
				Z2 = integrate(lattice, K2)
				@test _error(Z1, Z2, tol) < tol
				for i in 1:lattice.k-1
					for j in 1:lattice.k-1
						for f1 in (true, false), f2 in (true, false)
							for band in 1:lattice.bands
								g1 = gf(lattice_scaling, (i-1)*scaling+1, (j-1)*scaling+1, K1, f1=f1, f2=f2, band=band,Z=Z1)
								g2 = gf(lattice, i, j, K2, f1=f1, f2=f2, band=band,Z=Z2)
								@test _error(g1, g2, tol) < tol
							end		
						end
					end
				end
			end
		end
	end
end

@testset "build K-real time thermal state" begin
	tol = 1.0e-8
	U = 2.
	ϵ_d = 0.7
	δt = 0.02
	β = 1.
	D = 1.
	N = 3

	bath = fermionicbath(spectrum_func(D), β=β, μ=1)
	exact_model = SISB(bath, U=U, μ=ϵ_d)

	lattice = GrassmannLattice(δt=δt, N=N, bands=2, contour=:real)
	K = sysdynamics(lattice, exact_model)
	K = systhermalstate!(K, lattice, exact_model)
	for band in 1:lattice.bands
		K = boundarycondition(K, lattice, band=band)
	end
	n0 = occupation(lattice, K, band=1)

	for ordering in [A1A1a1a1B1B1b1b1(), A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2(), A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2()]
		lattice = GrassmannLattice(δt=δt, N=N, bands=2, contour=:real, ordering=ordering)
		K = accsysdynamics_fast(lattice, exact_model)
		K = systhermalstate!(K, lattice, exact_model)
		for band in 1:lattice.bands
			K = boundarycondition(K, lattice, band=band)
		end
		for band in 1:lattice.bands
			n2 = occupation(lattice, K, band=band)	
			@test norm(n2-n0) / norm(n0) < tol
		end		
	end

end

