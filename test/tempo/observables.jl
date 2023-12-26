println("------------------------------------")
println("|            Observables           |")
println("------------------------------------")




@testset "Occupation" begin
	rtol = 1.0e-7

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, order=1, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				B = randomgmps(Float64, length(lattice), D=6)
				cache = environments(lattice, A, B)
				Z = integrate(lattice, A, B)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = occupation(lattice, i, A, B, band=band, Z=Z)
						g2 = cached_occupation(lattice, i, A, B, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, order=2, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				B = randomgmps(Float64, length(lattice), D=6)
				cache = environments(lattice, A, B)
				Z = integrate(lattice, A, B)
				for band in 1:lattice.bands
					g1 = occupation(lattice, A, B, band=band, Z=Z)
					g2 = cached_occupation(lattice, A, B, cache=cache, band=band)
					@test abs(g1-g2)/abs(g1) < rtol
				end
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, order=1, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=4)
				AB = A * B
				cache = environments(lattice, A, B, C)
				Z = integrate(lattice, A, B, C)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = occupation(lattice, i, AB, C, band=band, Z=Z)
						g2 = cached_occupation(lattice, i, A, B, C, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, order=2, ordering=ordering)
				A = randomgmps(ComplexF64, length(lattice), D=2)
				B = randomgmps(ComplexF64, length(lattice), D=2)
				C = randomgmps(ComplexF64, length(lattice), D=4)
				AB = A * B
				cache = environments(lattice, A, B, C)
				Z = integrate(lattice, A, B, C)
				for band in 1:lattice.bands
					g1 = occupation(lattice, C, AB, band=band, Z=Z)
					g2 = cached_occupation(lattice, A, B, C, cache=cache, band=band)
					@test abs(g1-g2)/abs(g1) < rtol
				end
			end
		end
	end
end


@testset "Electric Current" begin
	rtol = 1.0e-7
	bath = fermionicbath(spectrum_func(1), β=1., μ=0)

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, order=1, ordering=ordering)
				corr = correlationfunction(bath, lattice)
				A = randomgmps(ComplexF64, length(lattice), D=4)
				B = randomgmps(ComplexF64, length(lattice), D=6)
				cache = environments(lattice, A, B)
				Z = integrate(lattice, A, B)
				for i in 2:lattice.k
					for band in 1:lattice.bands
						g1 = electriccurrent(lattice, corr, i, A, B, band=band, Z=Z)
						g2 = cached_electriccurrent(lattice, corr, i, A, B, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
						g3 = electriccurrent2(lattice, corr, i, A, B, band=band, Z=Z)
						@test abs(g1-g3)/abs(g1) < rtol
						g4 = cached_electriccurrent2(lattice, corr, i, A, B, cache=cache, band=band)
						@test abs(g1-g4)/abs(g1) < rtol
					end
				end
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, order=1, ordering=ordering)
				corr = correlationfunction(bath, lattice)
				A = randomgmps(ComplexF64, length(lattice), D=2)
				B = randomgmps(ComplexF64, length(lattice), D=2)
				C = randomgmps(ComplexF64, length(lattice), D=4)
				AB = A * B
				cache = environments(lattice, A, B, C)
				Z = integrate(lattice, A, B, C)
				for i in 2:lattice.k
					for band in 1:lattice.bands
						g1 = electriccurrent(lattice, corr, i, AB, C, band=band, Z=Z)
						g2 = cached_electriccurrent(lattice, corr, i, A, B, C, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
						g3 = electriccurrent2(lattice, corr, i, A, B, C, band=band, Z=Z)
						@test abs(g1-g3)/abs(g1) < rtol
						g4 = cached_electriccurrent2(lattice, corr, i, A, B, C, cache=cache, band=band)
						@test abs(g1-g4)/abs(g1) < rtol
					end
				end
			end
		end
	end

	bath = fermionicbath(spectrum_func(1), β=1., μ=1)

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, order=2, ordering=ordering)
				corr = correlationfunction(bath, lattice)
				A = randomgmps(ComplexF64, length(lattice), D=4)
				B = randomgmps(ComplexF64, length(lattice), D=6)
				cache = environments(lattice, A, B)
				Z = integrate(lattice, A, B)
				for band in 1:lattice.bands
					g1 = electriccurrent(lattice, corr, A, B, band=band, Z=Z)
					g2 = cached_electriccurrent(lattice, corr, A, B, cache=cache, band=band)
					@test abs(g1-g2)/abs(g1) < rtol
					g3 = electriccurrent2(lattice, corr, A, B, band=band, Z=Z)
					@test abs(g1-g3)/abs(g1) < rtol
					g4 = cached_electriccurrent2(lattice, corr, A, B, cache=cache, band=band)
					@test abs(g1-g4)/abs(g1) < rtol
				end
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, order=2, ordering=ordering)
				corr = correlationfunction(bath, lattice)
				A = randomgmps(ComplexF64, length(lattice), D=2)
				B = randomgmps(ComplexF64, length(lattice), D=2)
				C = randomgmps(ComplexF64, length(lattice), D=4)
				BC = B * C
				cache = environments(lattice, A, B, C)
				Z = integrate(lattice, A, B, C)
				for band in 1:lattice.bands
					g1 = electriccurrent(lattice, corr, A, BC, band=band, Z=Z)
					g2 = cached_electriccurrent(lattice, corr, A, B, C, cache=cache, band=band)
					@test abs(g1-g2)/abs(g1) < rtol
					g3 = electriccurrent2(lattice, corr, A, B, C, band=band, Z=Z)
					@test abs(g1-g3)/abs(g1) < rtol
					g4 = cached_electriccurrent2(lattice, corr, A, B, C, cache=cache, band=band)
					@test abs(g1-g4)/abs(g1) < rtol
				end
			end
		end
	end
end