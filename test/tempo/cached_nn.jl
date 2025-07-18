println("------------------------------------")
println("|             Cached nn            |")
println("------------------------------------")


@testset "Cached density-density correlation" begin

	rtol = 1.0e-5

	for N in 2:3
		for bands in 1:2
			for ordering in imag_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				cache = environments(lattice, A)
				Z = integrate(lattice, A)
				for i in 1:lattice.k-1, j in 1:lattice.k-1
					if j != i
						for band in 1:lattice.bands
							g1 = nn2(lattice, i, j, A, band=band, Z=Z)
							g2 = cached_nn2(lattice, i, j, A, cache=cache, band=band)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end				
			end
		end
	end


	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				cache = environments(lattice, A)
				Z = integrate(lattice, A)
				for i in 1:lattice.k-1, j in 1:lattice.k-1
					for b1 in branches(lattice), b2 in branches(lattice)
						if !((j == i) && (b1 == b2))
							g1 = nn2(lattice, i, j, A, b1=b1, b2=b2, Z=Z)
							g2 = cached_nn2(lattice, i, j, A, cache=cache, b1=b1, b2=b2)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end			
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in mixed_ac_grassmann_orderings
				lattice = GrassmannLattice(Nt=N, δt = 0.1, δτ=0.05, Nτ=2, bands=bands, contour=:mixed, ordering=ordering)
				A = randomgmps(scalartype(lattice), length(lattice), D=4)
				cache = environments(lattice, A)
				Z = integrate(lattice, A)
				for b1 in branches(lattice), b2 in branches(lattice)
					k1 = ifelse(b1==:τ, lattice.Nτ, lattice.Nt)
					k2 = ifelse(b2==:τ, lattice.Nτ, lattice.Nt)
					for i in 1:k1, j in 1:k2
						if !((j == i) && (b1 == b2))
							g1 = nn2(lattice, i, j, A, b1=b1, b2=b2, Z=Z)
							g2 = cached_nn2(lattice, i, j, A, cache=cache, b1=b1, b2=b2)
							@test abs(g1-g2)/abs(g1) < rtol							
						end
					end
				end

			end
		end
	end

end