println("------------------------------------")
println("|          Cached GF fast          |")
println("------------------------------------")

@testset "Cached GF fast" begin

	rtol = 1.0e-5

	for N in 2:3
		for bands in 1:2
			for ordering in imag_ac_grassmann_orderings
				# println("N=", N, " bands=", bands, " ordering=", ordering)
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				B = randomgmps(Float64, length(lattice), D=6)
				C = randomgmps(Float64, length(lattice), D=2)
				cache = environments(lattice, A)
				for band in 1:lattice.bands
					g1 = cached_Gτ(lattice, A, cache=cache, band=band)
					g2 = cached_Gτ_fast(lattice, A, cache=cache, band=band)
					@test norm(g1-g2)/norm(g1) < rtol
				end				
				cache = environments(lattice, A, B)
				for band in 1:lattice.bands
					g1 = cached_Gτ(lattice, A, B, cache=cache, band=band)
					g2 = cached_Gτ_fast(lattice, A, B, cache=cache, band=band)
					@test norm(g1-g2)/norm(g1) < rtol
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
				for band in 1:lattice.bands
					for f1 in (:+, ), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((f1 == f2) && (c1 == c2))
							# println("N=", N, " bands=", bands, " ordering=", ordering, " f1=", f1, " f2=", f2, " c1=", c1, " c2=", c2)
							g1 = [cached_Gt(lattice, i, 1, A, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2, band=band) for i in 1:lattice.k]
							g2 = cached_Gt_fast(lattice, A, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2, band=band)
							@test norm(g1-g2)/norm(g1) < rtol
						end
					end
				end

				B = randomgmps(Float64, length(lattice), D=6)
				cache = environments(lattice, A, B)
				
				for band in 1:lattice.bands
					for f1 in (:+, ), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
						if !((f1 == f2) && (c1 == c2))
							g1 = [cached_Gt(lattice, i, 1, A, B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2, band=band) for i in 1:lattice.k]
							g2 = cached_Gt_fast(lattice, A, B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2, band=band)
							@test norm(g1-g2)/norm(g1) < rtol
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

				for f1 in (:+, ), f2 in (:+, :-, :τ), c1 in (true, false), c2 in (true, false)
					if !((f1 == f2) && (c1 == c2))
						g1 = [cached_Gm(lattice, i, 1, A, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2) for i in 1:lattice.kt]
						g2 = cached_Gm_fast(lattice, A, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2) 
						@test norm(g1-g2)/norm(g1) < rtol
					end
				end

				B = randomgmps(scalartype(lattice), length(lattice), D=6)
				cache = environments(lattice, A, B)

				for f1 in (:τ,), f2 in (:+, :-, :τ,), c1 in (true, false), c2 in (true, false)
					if !((f1 == f2) && (c1 == c2))
						g1 = [cached_Gm(lattice, i, 1, A, B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2) for i in 1:lattice.kτ]
						g2 = cached_Gm_fast(lattice, A, B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2) 
						@test norm(g1-g2)/norm(g1) < rtol
					end
				end

			end
		end
	end

end
