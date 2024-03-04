println("------------------------------------")
println("|            BMPS GF               |")
println("------------------------------------")



@testset "BMPS GF" begin
	rtol = 1.0e-4

	for N in 2:3
		for bands in 1:2
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=2)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=4)
				D = randomgmps(Float64, length(lattice), D=2)

				Z1 = integrate(lattice, A, B, C, D)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = gf(lattice, i, A, B, C, D, band=band, Z=Z1)
						g2 = gf(lattice, i, A, B, C, D, band=band, Z=Z1, alg=BMPSIntegrate())
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end

				Z1 = integrate(lattice, [A, B, C, D], A, D)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = gf(lattice, i, [A, B, C, D], A, D, band=band, Z=Z1)
						g2 = gf(lattice, i, [A, B, C, D], A, D, band=band, Z=Z1, alg=BMPSIntegrate())
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end				
			end
		end
	end


	for N in 1:2
		for bands in 1:2
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, ordering=ordering)
				A = randomgmps(ComplexF64, length(lattice), D=4)
				B = randomgmps(ComplexF64, length(lattice), D=2)
				C = randomgmps(ComplexF64, length(lattice), D=2)
				D = randomgmps(ComplexF64, length(lattice), D=2)

				Z1 = integrate(lattice, A, B, C, D)
				for i in 1:lattice.k, j in 1:lattice.k
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = gf(lattice, i, j, A, B, C, D, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z1)
							g2 = gf(lattice, i, j, A, B, C, D, c1=c1, c2=c2, b1=f1, b2=f2, alg=BMPSIntegrate(), Z=Z1)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end

				Z1 = integrate(lattice, [A, B, C], D)
				for i in 1:lattice.k, j in 1:lattice.k
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = gf(lattice, i, j, [A, B, C], D, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z1)
							g2 = gf(lattice, i, j, [A, B, C], D, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z1, alg=BMPSIntegrate())
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end
			end
		end
	end

end