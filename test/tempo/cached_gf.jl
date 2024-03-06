println("------------------------------------")
println("|             Cached GF            |")
println("------------------------------------")

@testset "Cached GF" begin

	rtol = 1.0e-5

	for N in 2:3
		for bands in 1:2
			for ordering in imag_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				cache = environments(lattice, A)
				Z = integrate(lattice, A)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, A, band=band, Z=Z)
						g2 = cached_Gτ(lattice, i, A, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end				
				B = randomgmps(Float64, length(lattice), D=6)
				cache = environments(lattice, A, B)
				Z = integrate(lattice, A, B)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, A, B, band=band, Z=Z)
						g2 = cached_Gτ(lattice, i, A, B, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end
				C = randomgmps(Float64, length(lattice), D=2)
				cache = environments(lattice, [A, C], B)
				Z = integrate(lattice, [A, C], B)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, [A, C], B, band=band, Z=Z)
						g2 = cached_Gτ(lattice, i, [A,C], B, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
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
				for i in 1:lattice.k, j in 1:lattice.k
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = Gt(lattice, i, j, A, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z)
							g2 = cached_Gt(lattice, i, j, A, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end
				B = randomgmps(Float64, length(lattice), D=6)
				cache = environments(lattice, A, B)
				Z = integrate(lattice, A, B)
				for i in 1:lattice.k, j in 1:lattice.k
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = Gt(lattice, i, j, A, B, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z)
							g2 = cached_Gt(lattice, i, j, A, B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end
				C = randomgmps(Float64, length(lattice), D=2)
				cache = environments(lattice, [A, C], B)
				Z = integrate(lattice, [A, C], B)
				for i in 1:lattice.k, j in 1:lattice.k
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = Gt(lattice, i, j, [A,C], B, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z)
							g2 = cached_Gt(lattice, i, j, [A,C], B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2)
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
				for i in 1:lattice.kt, j in 1:lattice.kt
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = Gm(lattice, i, j, A, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z)
							g2 = cached_Gm(lattice, i, j, A, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end
				B = randomgmps(scalartype(lattice), length(lattice), D=6)
				cache = environments(lattice, A, B)
				Z = integrate(lattice, A, B)
				for i in 1:lattice.kτ, j in 1:lattice.kτ
					for f1 in (:τ,), f2 in (:τ,), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = Gm(lattice, i, j, A, B, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z)
							g2 = cached_Gm(lattice, i, j, A, B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end
				C = randomgmps(scalartype(lattice), length(lattice), D=2)
				cache = environments(lattice, [A, C], B)
				Z = integrate(lattice, [A, C], B)
				for i in 1:lattice.kt, j in 1:lattice.kt
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = Gm(lattice, i, j, [A,C], B, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z)
							g2 = cached_Gm(lattice, i, j, [A,C], B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end				
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in imag_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=4)
				cache = environments(lattice, A, B, C)
				Z = integrate(lattice, A, B, C)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, A, B, C, band=band, Z=Z)
						g2 = cached_Gτ(lattice, i, A, B, C, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end
			end
		end
	end

	for N in 2:3
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, ordering=ordering)
				A = randomgmps(ComplexF64, length(lattice), D=4)
				B = randomgmps(ComplexF64, length(lattice), D=2)
				C = randomgmps(ComplexF64, length(lattice), D=4)
				cache = environments(lattice, A, B, C)
				Z = integrate(lattice, A, B, C)
				for i in 1:lattice.k, j in 1:lattice.k
					for f1 in (:+, :-), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((i == j) && (f1 == f2) && (c1 == c2))
							g1 = Gt(lattice, i, j, A, B, C, c1=c1, c2=c2, b1=f1, b2=f2, Z=Z)
							g2 = cached_Gt(lattice, i, j, A, B, C, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2)
							@test abs(g1-g2)/abs(g1) < rtol
						end
					end
				end
			end
		end
	end
end

@testset "Cached GF: 4" begin
	rtol = 1.0e-4

	for N in 2:3
		for bands in 1:2
			for ordering in imag_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=2)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=2)
				D = randomgmps(Float64, length(lattice), D=4)

				AB = A * B
				CD = C * D 
				cache = environments(lattice, A, B, C, D)
				Z = integrate(lattice, AB, CD)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, AB, CD, band=band, Z=Z)
						g2 = cached_Gτ(lattice, i, A, B, C, D, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end
			end
		end
	end
end

@testset "Cached GF: 5" begin
	rtol = 1.0e-4

	for N in 2:3
		for bands in 1:2
			for ordering in imag_ac_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=2)
				B = randomgmps(Float64, length(lattice), D=2)
				C = randomgmps(Float64, length(lattice), D=2)
				D = randomgmps(Float64, length(lattice), D=2)
				E = randomgmps(Float64, length(lattice), D=2)

				ABC = A * B * C
				DE = D * E  
				cache = environments(lattice, E,D,C,B,A)
				Z = integrate(lattice, ABC, DE)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						g1 = Gτ(lattice, i, ABC, DE, band=band, Z=Z)
						g2 = cached_Gτ(lattice, i, E,D,C,B,A, cache=cache, band=band)
						@test abs(g1-g2)/abs(g1) < rtol
					end
				end
			end
		end
	end
end

# @testset "Cached GF: 6" begin
# 	rtol = 1.0e-4

# 	for N in 2:3
# 		for bands in 1:2
# 			for ordering in imag_ac_grassmann_orderings
# 				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
# 				A = randomgmps(Float64, length(lattice), D=2)
# 				B = randomgmps(Float64, length(lattice), D=2)
# 				C = randomgmps(Float64, length(lattice), D=2)
# 				D = randomgmps(Float64, length(lattice), D=2)
# 				E = randomgmps(Float64, length(lattice), D=2)
# 				F = randomgmps(Float64, length(lattice), D=2)

# 				ABC = A * B * C
# 				DEF = D * E * F 
# 				cache = environments(lattice, F,E,D,C,B,A)
# 				Z = integrate(lattice, ABC, DEF)
# 				for i in 1:lattice.k
# 					for band in 1:lattice.bands
# 						g1 = gf(lattice, i, ABC, DEF, band=band, Z=Z)
# 						g2 = cached_gf(lattice, i, F,E,D,C,B,A, cache=cache, band=band)
# 						@test abs(g1-g2)/abs(g1) < rtol
# 					end
# 				end
# 			end
# 		end
# 	end
# end

# @testset "Cached GF: 7" begin
# 	rtol = 1.0e-4

# 	for N in 2:3
# 		for bands in 1:2
# 			for ordering in imag_ac_grassmann_orderings
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
# 				cache = environments(lattice, A, B, C, D, E, F, G)
# 				Z = integrate(lattice, ABC, DEFG)
# 				for i in 1:lattice.k
# 					for band in 1:lattice.bands
# 						g1 = gf(lattice, i, ABC, DEFG, band=band, Z=Z)
# 						g2 = cached_gf(lattice, i, A, B, C, D, E, F, G, cache=cache, band=band)
# 						@test abs(g1-g2)/abs(g1) < rtol
# 					end
# 				end
# 			end
# 		end
# 	end
# end