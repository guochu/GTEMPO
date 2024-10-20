println("------------------------------------")
println("|    PartialIF-Retarded Interact   |")
println("------------------------------------")


@testset "InfluenceFunctional: imaginary time" begin
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)
	ordering = A1A1B1B1()
	for bands in (1,2,3)
		for N in (2,3)

			lattice = GrassmannLattice(N=N, δτ=0.05, bands=bands, contour=:imag, order=1, ordering=ordering)
			η = randn(Float64, lattice.k-1, lattice.k-1)
			for i in 1:lattice.k-1
				for band in 1:lattice.bands
					mps1 = GrassmannMPS(length(lattice))
					for j in 1:lattice.k-1
						if i == j
							pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
							t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
						else
							pos1a, pos1b = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
							pos2a, pos2b = index(lattice, j+1, conj=true, band=band), index(lattice, j, conj=false, band=band)
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