println("------------------------------------")
println("|          Cached GF fast2         |")
println("------------------------------------")

@testset "Cached GF fast2" begin
	rtol = 1.0e-5

	for N in 3:4
		for bands in 1:2
			for ordering in real_ac_grassmann_orderings
				for idx0 in 1;2

					lattice = GrassmannLattice(N=N, δt = 0.1, bands=bands, contour=:real, ordering=ordering)
					A = randomgmps(Float64, length(lattice), D=4)
					cache = environments2(lattice, A)
					for band in 1:lattice.bands
						for f1 in (:+, ), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
							if !((f1 == f2) && (c1 == c2))
								g1 = [cached_Gt(lattice, i, idx0, A, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2, band=band) for i in idx0:lattice.k]
								g2 = cached_Gt_fast(idx0, lattice, A, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2, band=band)
								@test norm(g1-g2)/norm(g1) < rtol
							end
						end

						g1 = [cached_greater(lattice, i, idx0, A, cache=cache, band=band) for i in idx0:lattice.k]
						g2 = cached_greater_fast(idx0, lattice, A, cache=cache, band=band)
						@test norm(g1-g2)/norm(g1) < rtol
						g1 = [cached_lesser(lattice, idx0, i, A, cache=cache, band=band) for i in idx0:lattice.k]
						g2 = cached_lesser_fast(idx0, lattice, A, cache=cache, band=band)
						@test norm(g1-g2)/norm(g1) < rtol
						
					end

					B = randomgmps(Float64, length(lattice), D=6)
					cache = environments2(lattice, A, B)
					
					for band in 1:lattice.bands
						for f1 in (:+, ), f2 in (:+, :-), c1 in (true, false), c2 in (true, false)
							if !((f1 == f2) && (c1 == c2))
								g1 = [cached_Gt(lattice, i, idx0, A, B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2, band=band) for i in idx0:lattice.k]
								g2 = cached_Gt_fast(idx0, lattice, A, B, cache=cache, c1=c1, c2=c2, b1=f1, b2=f2, band=band)
								@test norm(g1-g2)/norm(g1) < rtol
							end
						end

						g1 = [cached_greater(lattice, i, idx0, A, B, cache=cache, band=band) for i in idx0:lattice.k]
						g2 = cached_greater_fast(idx0, lattice, A, B, cache=cache, band=band)
						@test norm(g1-g2)/norm(g1) < rtol
						g1 = [cached_lesser(lattice, idx0, i, A, B, cache=cache, band=band) for i in idx0:lattice.k]
						g2 = cached_lesser_fast(idx0, lattice, A, B, cache=cache, band=band)
						@test norm(g1-g2)/norm(g1) < rtol	
					end
				end
			end
		end
	end

end
