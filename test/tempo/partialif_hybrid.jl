println("------------------------------------")
println("|         PartialIF-Hybrid         |")
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

						mps3 = partialif_hybrid(lattice, i, η[i, :], band=band)	
						@test bond_dimension(mps3) == 2
						@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5

						mps4 = partialif_hybrid_naive(lattice, i, η[i, :], band=band, trunc=trunc)	
						@test distance(mps4, mps3) / norm(mps3) <= 1.0e-5
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
	f = spectrum(J)
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)
	for bands in (1, 2)
		for N in (1, 2, 3)
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt=0.05, bands=bands, contour=:real, order=1, ordering=ordering)
				corr = fermionic_Δt(f, β=β, N=lattice.N, t=lattice.t)
				for fi in (:+, :-)
					for fj in (:+, :-)
						η = branch(corr, fi, fj)
						@test size(η, 1) == lattice.k
						for i in 1:lattice.k
							for band in 1:lattice.bands
								mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
								for j in 1:lattice.k
									pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=fj, band=band)
									t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
									mps1 = t * mps1
									canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))
								end

								mps3 = partialif_hybrid(lattice, i, view(η, i, 1:lattice.k), band=band, b1=fi, b2=fj)
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
				corr = fermionic_Δt(f, β=β, N=lattice.N, t=lattice.t)
				for fi in (:+, :-)
					η₁, η₂ = branch(corr, fi, :+), branch(corr, fi, :-)
					for i in 1:lattice.k
						for band in 1:lattice.bands
							mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
							for j in 1:lattice.k
								pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=:+, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η₁[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))			
								pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=:-, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η₂[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))															
							end

							mps3 = partialif_hybrid(lattice, i, view(η₁, i, 1:lattice.k), view(η₂, i, 1:lattice.k), band=band, b1=fi)
							@test bond_dimension(mps3) == 2
							@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5

							mps4 = partialif_hybrid_naive(lattice, i, view(η₁, i, 1:lattice.k), view(η₂, i, 1:lattice.k), band=band, b1=fi, trunc=trunc)
							@test distance(mps4, mps3) / norm(mps3) <= 1.0e-5
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
				for fi in (:+, :-), fj in (:+, :-)

					for j in 1:lattice.k
						for band in 1:lattice.bands
							mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
							for i in 1:lattice.k
								pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=fj, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))		
							end

							mps3 = partialif_hybrid(lattice, η[1:lattice.k, j], j, band=band, b1=fi, b2=fj)
							@test bond_dimension(mps3) == 2
							@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
						end
					end
				end
				η₂ = randn(ComplexF64, lattice.k, lattice.k)
				for fi in (:+, :-)
					for i in 1:lattice.k
						for band in 1:lattice.bands
							mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
							for j in 1:lattice.k
								pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=:+, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))
								pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=:-, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η₂[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))
							end
							mps2 = partialif_hybrid(lattice, i, η[i, 1:lattice.k], η₂[i, 1:lattice.k], band=band, b1=fi)
							@test distance(mps1, mps2) / norm(mps1) <= 1.0e-5

							mps3 = partialif_hybrid_naive(lattice, i, η[i, 1:lattice.k], η₂[i, 1:lattice.k], band=band, b1=fi, trunc=trunc)
							@test distance(mps3, mps2) / norm(mps2) <= 1.0e-5
						end
					end
				end
				for fj in (:+, :-)
					for j in 1:lattice.k
						for band in 1:lattice.bands
							mps1 = GrassmannMPS(scalartype(lattice), length(lattice))
							for i in 1:lattice.k
								pos1, pos2 = index(lattice, i, conj=true, branch=:+, band=band), index(lattice, j, conj=false, branch=fj, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))		
								pos1, pos2 = index(lattice, i, conj=true, branch=:-, band=band), index(lattice, j, conj=false, branch=fj, band=band)
								t = exp(GTerm(pos1, pos2, coeff=η₂[i, j]))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))		
							end
							mps3 = partialif_hybrid(lattice, η[1:lattice.k, j], η₂[1:lattice.k, j], j, band=band, b2=fj)
							@test bond_dimension(mps3) == 2
							@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
						end
					end
				end
			end
		end
	end

end



@testset "InfluenceFunctional: mixed time" begin
	D = 1
	β=0.2
	Nτ = 2
	δτ = β / Nτ
	J(ω) = (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1
	f = spectrum(J, lb = -D, ub = D)
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

	for bands in (1, 2)
		for N in (1, 2, 3)
			for ordering in mixed_grassmann_orderings
				lattice = GrassmannLattice(Nt=N, δt=0.05, δτ=δτ, Nτ=Nτ, bands=bands, contour=:mixed, order=1, ordering=ordering)
				corr = fermionic_Δm(f, β=β, Nt=lattice.Nt, t=lattice.t, Nτ=lattice.Nτ)
				for fi in (:+, :-)
					for i in 1:lattice.kt
						cols_f = [index(corr, i, j, b1=fi, b2=:+) for j in 1:lattice.kt]
						cols_b = [index(corr, i, j, b1=fi, b2=:-) for j in 1:lattice.kt]
						cols_i = [index(corr, i, j, b1=fi, b2=:τ) for j in 1:lattice.Nτ]

						for band in 1:lattice.bands
							mps1 = vacuumstate(lattice)
							for j in 1:lattice.kt
								pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=:+, band=band)
								t = exp(GTerm(pos1, pos2, coeff=index(corr, i, j, b1=fi, b2=:+)))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))			
								pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=:-, band=band)
								t = exp(GTerm(pos1, pos2, coeff=index(corr, i, j, b1=fi, b2=:-)))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))															
							end
							for j in 1:lattice.Nτ
								pos1, pos2 = index(lattice, i, conj=true, branch=fi, band=band), index(lattice, j+1, conj=false, branch=:τ, band=band)
								t = exp(GTerm(pos1, pos2, coeff=index(corr, i, j, b1=fi, b2=:τ)))
								mps1 = t * mps1
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))			
							end

							mps3 = partialif_hybrid(lattice, i, cols_f, cols_b, cols_i, band=band, b1=fi)
							@test bond_dimension(mps3) == 2
							@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5

							mps4 = partialif_hybrid_naive(lattice, i, cols_f, cols_b, cols_i, band=band, b1=fi, trunc=trunc)
							@test distance(mps4, mps3) / norm(mps3) <= 1.0e-5
						end
					end						
				end
				fi = :τ
				for i in 1:lattice.Nτ
					cols_f = [index(corr, i, j, b1=fi, b2=:+) for j in 1:lattice.kt]
					cols_b = [index(corr, i, j, b1=fi, b2=:-) for j in 1:lattice.kt]
					cols_i = [index(corr, i, j, b1=fi, b2=:τ) for j in 1:lattice.Nτ]

					for band in 1:lattice.bands
						mps1 = vacuumstate(lattice)
						for j in 1:lattice.kt
							pos1, pos2 = index(lattice, i+1, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=:+, band=band)
							t = exp(GTerm(pos1, pos2, coeff=index(corr, i, j, b1=fi, b2=:+)))
							mps1 = t * mps1
							canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))			
							pos1, pos2 = index(lattice, i+1, conj=true, branch=fi, band=band), index(lattice, j, conj=false, branch=:-, band=band)
							t = exp(GTerm(pos1, pos2, coeff=index(corr, i, j, b1=fi, b2=:-)))
							mps1 = t * mps1
							canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))															
						end
						for j in 1:lattice.Nτ
							pos1, pos2 = index(lattice, i+1, conj=true, branch=fi, band=band), index(lattice, j+1, conj=false, branch=:τ, band=band)
							t = exp(GTerm(pos1, pos2, coeff=index(corr, i, j, b1=fi, b2=:τ)))
							mps1 = t * mps1
							canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))			
						end

						mps3 = partialif_hybrid(lattice, i, cols_f, cols_b, cols_i, band=band, b1=fi)
						@test bond_dimension(mps3) == 2
						@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5

						mps4 = partialif_hybrid_naive(lattice, i, cols_f, cols_b, cols_i, band=band, b1=fi, trunc=trunc)
						@test distance(mps4, mps3) / norm(mps3) <= 1.0e-5
					end
				end						
			end
		end
	end

end


