# println("------------------------------------")
# println("|         PartialIF-Retarded       |")
# println("------------------------------------")


# @testset "InfluenceFunctional for Retarded Interact: imaginary time" begin
# 	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

# 	# one band
# 	for N in (2,3)
# 		lattice = FockLattice(N=N, δτ=0.05, contour=:imag, order=1)
# 		η = randn(Float64, lattice.N, lattice.N)
# 		band = 1
# 		for i in 1:lattice.N
# 			println("N = ", N, ", i = ", i)
# 			pos1 = index(lattice, i, band=band)
# 			mps1 = vacuumstate(lattice)
# 			coefs = scalartype(lattice)[]
# 			cols = Int[]
# 			for j in 1:lattice.N
# 				pos2 = index(lattice, j, band=band)
# 				coef = exp(η[i, j]) - 1
# 				push!(cols, pos2)
# 				push!(coefs, coef)
# 				if pos1 == pos2
# 					t = exp(NTerm(pos1, coeff=coef))
# 				else
# 					t = exp(NTerm(pos1, pos2, coeff=coef))
# 				end
# 				apply!(t, mps1)
# 				canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))
# 			end	
# 			mps2 = partialdensempo(pos1, cols, coefs) * vacuumstate(lattice)

# 			@test bond_dimension(mps2) == 2
# 			@test distance(mps1, mps2) / norm(mps1) <= 1.0e-5

# 		end
# 	end		
# end
