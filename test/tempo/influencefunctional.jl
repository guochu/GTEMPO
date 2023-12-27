println("------------------------------------")
println("|       InfluenceFunctional        |")
println("------------------------------------")


@testset "InfluenceFunctional: imaginary time" begin
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)
	for bands in (1,2,3)
		for N in (1,2,3)
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.05, bands=bands, contour=:imag, order=1, ordering=ordering)
				η = randn(Float64, lattice.k, lattice.k)
				for i in 1:lattice.k
					for band in 1:lattice.bands
						mps1 = GrassmannMPS(length(lattice))
						for j in 1:lattice.k
							pos1, pos2 = index(lattice, i, conj=true, band=band), index(lattice, j, conj=false, band=band)
							t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
							mps1 = t * mps1
							canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))
						end	

						mps3 = partialinfluencefunctional(lattice, i, η[i, :], band=band)	
						@test bond_dimension(mps3) == 2
						@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
					end
				end
			end
		end		
	end
end

@testset "InfluenceFunctional: real time" begin
	D = 1
	β=0.2
	J(ε) = D/(ε^2+D^2)/pi
	f = SpectrumFunction(J)
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)
	for bands in (1, 2)
		for N in (1, 2, 3)
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt=0.05, bands=bands, contour=:real, order=1, ordering=ordering)
				corr = Gt(f, β=β, N=lattice.N, t=lattice.t)
				for fi in (true, false)
					for fj in (true, false)
						η = branch(corr, fi, fj)
						@test size(η, 1) == lattice.k
						for i in 1:lattice.k
							for band in 1:lattice.bands
								mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
								for j in 1:lattice.k
									pos1, pos2 = index(lattice, i, conj=true, forward=fi, band=band), index(lattice, j, conj=false, forward=fj, band=band)
									t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
									mps1 = t * mps1
									canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))
								end

								mps3 = partialinfluencefunctional(lattice, i, view(η, i, 1:lattice.k), band=band, fi=fi, fj=fj)
								@test bond_dimension(mps3) == 2
								@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
							end
						end						
					end
				end
			end
		end
	end

	for bands in (1, 2)
		for N in (1, 2, 3)
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt=0.05, bands=bands, contour=:real, order=1, ordering=ordering)
				corr = Gt(f, β=β, N=lattice.N, t=lattice.t)
				for fi in (true, false)
					η₁, η₂ = branch(corr, fi, true), branch(corr, fi, false)
					for i in 1:lattice.k
						for band in 1:lattice.bands
							mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
							for j in 1:lattice.k
								pos1, pos2 = index(lattice, i, conj=true, forward=fi, band=band), index(lattice, j, conj=false, forward=true, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η₁[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))			
								pos1, pos2 = index(lattice, i, conj=true, forward=fi, band=band), index(lattice, j, conj=false, forward=false, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η₂[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))															
							end

							mps3 = partialinfluencefunctional(lattice, i, view(η₁, i, 1:lattice.k), view(η₂, i, 1:lattice.k), band=band, fi=fi)
							@test bond_dimension(mps3) == 2
							@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
						end
					end						
				end
			end
		end
	end

	for bands in (1, 2)
		for N in (1, 2, 3)
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt=0.05, bands=bands, contour=:real, order=1, ordering=ordering)
				η = randn(ComplexF64, lattice.k, lattice.k)
				for fi in (true, false), fj in (true, false)

					for j in 1:lattice.k
						for band in 1:lattice.bands
							mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
							for i in 1:lattice.k
								pos1, pos2 = index(lattice, i, conj=true, forward=fi, band=band), index(lattice, j, conj=false, forward=fj, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))		
							end

							mps3 = partialinfluencefunctional(lattice, η[1:lattice.k, j], j, band=band, fi=fi, fj=fj)
							@test bond_dimension(mps3) == 2
							@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
						end
					end
				end
				η₂ = randn(ComplexF64, lattice.k, lattice.k)
				for fi in (true, false)
					for i in 1:lattice.k
						for band in 1:lattice.bands
							mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
							for j in 1:lattice.k
								pos1, pos2 = index(lattice, i, conj=true, forward=fi, band=band), index(lattice, j, conj=false, forward=true, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))
								pos1, pos2 = index(lattice, i, conj=true, forward=fi, band=band), index(lattice, j, conj=false, forward=false, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η₂[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))
							end
							mps2 = partialinfluencefunctional(lattice, i, η[i, 1:lattice.k], η₂[i, 1:lattice.k], band=band, fi=fi)
							@test distance(mps1, mps2) / norm(mps1) <= 1.0e-5
						end
					end
				end
				for fj in (true, false)
					for j in 1:lattice.k
						for band in 1:lattice.bands
							mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
							for i in 1:lattice.k
								pos1, pos2 = index(lattice, i, conj=true, forward=true, band=band), index(lattice, j, conj=false, forward=fj, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))		
								pos1, pos2 = index(lattice, i, conj=true, forward=false, band=band), index(lattice, j, conj=false, forward=fj, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η₂[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))		
							end
							mps3 = partialinfluencefunctional(lattice, η[1:lattice.k, j], η₂[1:lattice.k, j], j, band=band, fj=fj)
							@test bond_dimension(mps3) == 2
							@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
						end
					end
				end
			end
		end
	end

end

