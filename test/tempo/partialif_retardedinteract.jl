println("------------------------------------")
println("|    PartialIF-Retarded Interact   |")
println("------------------------------------")


@testset "Retarded Interact IF: imaginary time" begin
	trunc = truncdimcutoff(D=300, ϵ=1.0e-8, add_back=0)
	for bands in (1,2,3)
		for N in (2,3)
			for ordering in imag_grassmann_orderings

				# println("bands = ", bands, ", N = ", N)
				lattice = GrassmannLattice(N=N, δτ=0.05, bands=bands, contour=:imag, order=1, ordering=ordering)
				η = randn(Float64, lattice.k-1, lattice.k-1)
				for i in 1:lattice.k-1, band in 1:lattice.bands
					mps1 = vacuumstate(lattice)
					for j in 1:lattice.k-1, band2 in 1:lattice.bands
						if (i == j) && (band == band2)
							pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
							t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
						else
							pos1a, pos1b = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
							pos2a, pos2b = index(lattice, j+1, conj=true, band=band2), index(lattice, j, conj=false, band=band2)
							t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=η[i, j]))
						end
						mps1 = t * mps1
						canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))
					end	

					mps3 = partialif_retardedinteract(lattice, i, η[i, :], band=band)	
					# @test bond_dimension(mps3) == 2
					@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
				end

			end

		end		
	end
end

@testset "Retarded Interact IF: real time" begin
	D = 1
	β=0.2
	J(ε) = D/(ε^2+D^2)/pi
	f = SpectrumFunction(J)
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

	for bands in (1, 2)
		for N in (2, 3, 4)
			for ordering in real_grassmann_orderings
				# println("bands=", bands, ", N=", N)
				lattice = GrassmannLattice(N=N, δt=0.05, bands=bands, contour=:real, order=1, ordering=ordering)
				corr = fermionic_Ct(f, β=β, N=lattice.N, t=lattice.t)
				k = lattice.k-1
				for fi in (:+, :-)
					η₁, η₂ = branch(corr, fi, :+), branch(corr, fi, :-)
					for i in 1:k, band in 1:lattice.bands
						mps1 = vacuumstate(lattice)
						for j in 1:k, band2 in 1:lattice.bands
							pos1a, pos1b = index(lattice, i+1, conj=true, band=band, branch=fi), index(lattice, i, conj=false, band=band, branch=fi)
							for fj in (:+, :-)
								η = (fj == :+) ? η₁ : η₂
								if (i == j) && (band == band2) && (fi == fj)
									t = exp(GTerm(pos1a, pos1b, coeff=η[i, j]))
								else
									pos2a, pos2b = index(lattice, j+1, conj=true, band=band2, branch=fj), index(lattice, j, conj=false, band=band2, branch=fj)
									t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=η[i, j]))
								end
								apply!(t, mps1)
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))	
							end
						end

						mps3 = partialif_retardedinteract(lattice, i, view(η₁, i, 1:k), view(η₂, i, 1:k), band=band, b1=fi)
						# println("bond dimension ", bond_dimension(mps1), " ", bond_dimension(mps3)) 
						@test distance(mps1, mps3) / norm(mps1) <= 1.0e-5
						# println("i=", i, " error=", distance(mps1, mps3) / norm(mps1))
					end						
				end
			end
		end
	end

end
