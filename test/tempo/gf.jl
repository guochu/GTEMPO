println("------------------------------------")
println("|                GF                |")
println("------------------------------------")


@testset "GF" begin
	rtol = 1.0e-7

	for N in 2:3
		for bands in 1:2
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=2)
				AB = A * B
				ABC = A * B * C
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, AB, band=band)
						g2 = Gτ(lattice, i, A, B, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
						g3 = Gτ(lattice, i, B, A, band=band)
						@test abs(g1-g3)/abs(g1) < rtol

						g1 = Gτ(lattice, i, ABC, band=band)
						g2 = Gτ(lattice, i, A, B, C, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
						g3 = Gτ(lattice, i, C, A, B, band=band)
						@test abs(g1-g3)/abs(g1) < rtol
						g4 = Gτ(lattice, i, C, AB, band=band)
						@test abs(g1-g4)/abs(g1) < rtol


					end
				end
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=2)
				AB = A + B
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, AB, band=band)
						g2 = Gτ(lattice, i, [A, B], band=band)
						@test abs(g1-g2)/abs(g1) < rtol
						g3 = Gτ(lattice, i, [B, A], band=band)
						@test abs(g1-g3)/abs(g1) < rtol

						g1 = Gτ(lattice, i, AB, C, band=band)
						g2 = Gτ(lattice, i, [A, B], C, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
						g3 = Gτ(lattice, i, [B, A], C, band=band)
						@test abs(g1-g3)/abs(g1) < rtol
					end
				end
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, ordering=ordering)
				A = randomgmps(ComplexF64, length(lattice), D=4)
				B = randomgmps(ComplexF64, length(lattice), D=2)
				C = randomgmps(ComplexF64, length(lattice), D=4)
				AB = A + B
				Z1 = integrate(lattice, AB, C)
				Z2 = integrate(lattice, [A, B], C)
				@test abs(Z1-Z2)/abs(Z1) < rtol
				for i in 1:lattice.k, j in 1:lattice.k
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = Gt(lattice, i, j, AB, C, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z1)
							g2 = Gt(lattice, i, j, [A, B], C, Z=Z2, c1=c1, c2=c2, b1=f1, b2=f2)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end
			end
		end
	end


end

@testset "GF: 4" begin
	rtol = 1.0e-4

	for N in 2:3
		for bands in 1:2
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=2)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=4)
				D = randomgmps(Float64, length(lattice), D=2)

				AB = A * B
				CD = C * D 
				Z1 = integrate(lattice, AB, CD)
				Z2 = integrate(lattice, A, B, C, D)
				@test abs(Z1-Z2)/abs(Z1) < rtol
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, AB, CD, band=band, Z=Z1)
						g2 = Gτ(lattice, i, A, B, C, D, band=band, Z=Z2)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end
			end
		end
	end
end

@testset "GF: 5" begin
	rtol = 1.0e-4

	for N in 2:3
		for bands in 1:2
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=2)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=2)
				D = randomgmps(Float64, length(lattice), D=2)
				E = randomgmps(Float64, length(lattice), D=2)

				ABC = A * B * C
				DE = D * E 
				Z1 = integrate(lattice, ABC, DE)
				Z2 = integrate(lattice, A, B, E, D, C)
				@test abs(Z1-Z2)/abs(Z1) < rtol
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, ABC, DE, band=band, Z=Z1)
						g2 = Gτ(lattice, i, A, B, E, D, C, band=band, Z=Z2)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end
			end
		end
	end
end

# no longer supported
# @testset "GF: 6" begin
# 	rtol = 1.0e-4

# 	for N in 2:3
# 		for bands in 1:2
# 			for ordering in imag_grassmann_orderings
# 				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
# 				A = randomgmps(Float64, length(lattice), D=2)
# 				B = randomgmps(Float64, length(lattice), D=2)
# 				C = randomgmps(Float64, length(lattice), D=2)
# 				D = randomgmps(Float64, length(lattice), D=2)
# 				E = randomgmps(Float64, length(lattice), D=2)
# 				F = randomgmps(Float64, length(lattice), D=2)

# 				ABC = A * B * C
# 				DEF = D * E * F 
# 				Z1 = integrate(lattice, ABC, DEF)
# 				Z2 = integrate(lattice, A, B, F, E, D, C)
# 				@test abs(Z1-Z2)/abs(Z1) < rtol
# 				for i in 1:lattice.k
# 					for band in 1:lattice.bands
# 						g1 = gf(lattice, i, ABC, DEF, band=band, Z=Z1)
# 						g2 = gf(lattice, i, A, B, F, E, D, C, band=band, Z=Z2)
# 						@test abs(g1-g2)/abs(g1) < rtol
# 					end
# 				end
# 			end
# 		end
# 	end
# end

# @testset "GF: 7" begin
# 	rtol = 1.0e-4

# 	for N in 2:3
# 		for bands in 1:2
# 			for ordering in imag_grassmann_orderings
# 				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
# 				A = randomgmps(Float64, length(lattice), D=2)
# 				B = randomgmps(Float64, length(lattice), D=2)
# 				C = randomgmps(Float64, length(lattice), D=2)
# 				D = randomgmps(Float64, length(lattice), D=2)
# 				E = randomgmps(Float64, length(lattice), D=2)
# 				F = randomgmps(Float64, length(lattice), D=2)
# 				G = randomgmps(Float64, length(lattice), D=2)

# 				ABC = A * B * C
# 				DEFG = D * E * F * G
# 				Z1 = integrate(lattice, ABC, DEFG)
# 				Z2 = integrate(lattice, A, B, C, D, E, F, G)
# 				@test abs(Z1-Z2)/abs(Z1) < rtol
# 				for i in 1:lattice.k
# 					for band in 1:lattice.bands
# 						g1 = gf(lattice, i, ABC, DEFG, band=band, Z=Z1)
# 						g2 = gf(lattice, i, A, B, C, D, E, F, G, band=band, Z=Z2)
# 						@test abs(g1-g2)/abs(g1) < rtol
# 					end
# 				end
# 			end
# 		end
# 	end
# end
