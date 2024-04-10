println("------------------------------------")
println("|            GrassmannMPS          |")
println("------------------------------------")

@testset "GrassmannMPS: arithmetic and canonicalize" begin
	L = 6
	D = 6
	for T in (Float64, ComplexF64)
		psi = randomgmps(T, L, D=D)
		@test scalartype(psi) == T
		@test space_l(psi) == oneunit(grassmannpspace())
		@test space_r(psi) == oneunit(grassmannpspace())'

		@test bond_dimension(psi) <= D
		psi1 = leftorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = 1.0e-7
		@test distance(psi, psi1) < 1.0e-7

		psi1 = rightorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = 1.0e-7
		@test distance(psi, psi1) < 1.0e-7

		psi1 = leftorth!(deepcopy(psi), alg = Orthogonalize(QR(), normalize=true))
		@test isleftcanonical(psi1)
		psi1 = rightorth!(deepcopy(psi), alg = Orthogonalize(SVD(), normalize=true))
		@test isrightcanonical(psi1)
		psi1 = canonicalize!(deepcopy(psi), alg = Orthogonalize(SVD(), normalize=true))
		@test iscanonical(psi1)
		@test norm(2 * psi1) ≈ 2
		@test norm(psi1 / 2) ≈ 0.5
		@test norm(psi1 - psi1) ≈ 0. atol = 1.0e-7
		@test distance(psi, psi) ≈ 0. atol = 1.0e-7

		psi1 = canonicalize!(deepcopy(psi), alg=Orthogonalize(trunc=NoTruncation(), normalize=false))
		@test norm(psi) ≈ norm(psi1) atol = 1.0e-7
		@test distance(psi, psi1) < 1.0e-7

		# psi1 = randomgmps(T, L, D=2)
		# psi2 = randomgmps(T, L, D=4)

		# psi3 = psi1 * psi2
		# _n = norm(psi3)

		# canonicalize!(psi1)
		# canonicalize!(psi2)
		# psi5 = psi1 * psi2
		# @test distance(psi3, psi5) / _n < 1.0e-7

		# psi = randomgmps(T, L, D=D)
		# psi2 = increase_bond!(psi, 10)

		# @test distance(psi, psi2) / norm(psi) < 1.0e-7
	end

	


	
end

@testset "GrassmannMPS: multiplications" begin
	L = 6
	chi = 20
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	alg1 = SVDCompression(trunc)
	alg2 = DMRGMult1(trunc, initguess=:svd)
	alg3 = DMRGMult1(trunc, initguess=:rand, maxiter=10)
	alg4 = DMRGMult1(trunc, initguess=:pre, maxiter=10)
	alg5 = DMRGMult2(trunc, initguess=:svd)
	alg6 = DMRGMult2(trunc, initguess=:rand, maxiter=10)
	alg7 = DMRGMult2(trunc, initguess=:pre, maxiter=10)
	algs = [alg1, alg2, alg3, alg4, alg5, alg6, alg7]
	tol = 1.0e-7
	for T in (Float64, ComplexF64)
		psi1 = randomgmps(T, L, D=4)
		psi2 = randomgmps(T, L, D=4)

		psi3 = psi1 * psi2
		_n = norm(psi3)
		for alg in algs
			psi4 = mult(psi1, psi2, alg)
			@test distance(psi3, psi4) / _n < tol
		end

		canonicalize!(psi1)
		canonicalize!(psi2)
		psi5 = psi1 * psi2
		@test distance(psi3, psi5) / _n < tol

		for alg in algs
			psi4 = mult(psi1, psi2, alg)
			@test distance(psi3, psi4) / _n < tol
		end
	end
end


@testset "GrassmannMPS: ordering conversion" begin
	for N in (1, 2,3)
		for bands in (1,2,3)
			lattice = GrassmannLattice(δτ=0.1, N=N, bands=bands, contour=:imag, ordering=ABBA())
			K1 = randomgmps(scalartype(lattice), length(lattice), D=6)
			lattice2, K2 = toadjacentordering(lattice, K1)
			Z1 = integrate(lattice, K1)
			Z2 = integrate(lattice2, K2)
			@test abs((Z1-Z2) / Z1) <= 1.0e-5

			for ordering in (A1a1B1b1b1B1a1A1(), A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2())
				lattice = GrassmannLattice(δt=0.1, N=N, bands=bands, contour=:real, ordering=ordering)
				K1 = randomgmps(scalartype(lattice), length(lattice), D=6)
				lattice2, K2 = toadjacentordering(lattice, K1)
				Z1 = integrate(lattice, K1)
				Z2 = integrate(lattice2, K2)
				@test abs((Z1-Z2) / Z1) <= 1.0e-5
			end
		end
	end
end

@testset "GrassmannMPS: addition and integration" begin
	tol = 1.0e-8
	for N in (1, 2)
		for bands in (1,2,3)
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(scalartype(lattice), length(lattice), D=4)
				B = randomgmps(scalartype(lattice), length(lattice), D=6)
				C = randomgmps(scalartype(lattice), length(lattice), D=2)
				D = randomgmps(scalartype(lattice), length(lattice), D=2)

				Z1 = integrate(lattice, A + B)
				Z2 = integrate(lattice, [A, B])
				@test abs((Z1-Z2) / Z1) <= tol
				Z1 = integrate(lattice, A + B + C)
				Z2 = integrate(lattice, [A, B, C])
				@test abs((Z1-Z2) / Z1) <= tol
				Z1 = integrate(lattice, A + B, C)
				Z2 = integrate(lattice, [A, B], C)
				@test abs((Z1-Z2) / Z1) <= tol
				Z1 = integrate(lattice, A + B, C, D)
				Z2 = integrate(lattice, [A, B], C, D)
				@test abs((Z1-Z2) / Z1) <= tol
			end
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt=0.1, bands=bands, contour=:real, ordering=ordering)
				A = randomgmps(scalartype(lattice), length(lattice), D=4)
				B = randomgmps(scalartype(lattice), length(lattice), D=6)
				C = randomgmps(scalartype(lattice), length(lattice), D=2)
				D = randomgmps(scalartype(lattice), length(lattice), D=2)

				Z1 = integrate(lattice, A + B)
				Z2 = integrate(lattice, [A, B])
				@test abs((Z1-Z2) / Z1) <= tol
				Z1 = integrate(lattice, A + B + C)
				Z2 = integrate(lattice, [A, B, C])
				@test abs((Z1-Z2) / Z1) <= tol
				Z1 = integrate(lattice, A + B, C)
				Z2 = integrate(lattice, [A, B], C)
				@test abs((Z1-Z2) / Z1) <= tol
				Z1 = integrate(lattice, A + B, C, D)
				Z2 = integrate(lattice, [A, B], C, D)
				@test abs((Z1-Z2) / Z1) <= tol
			end
		end		
	end

end

@testset "GrassmannMPS: multiplication and integration" begin
	tol = 1.0e-5
	for N in (1, 2,3)
		for bands in (1,2,3)
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=4)
				B = randomgmps(Float64, length(lattice), D=6)
				C = A * B
				Z1 = integrate(lattice, C)
				Z2 = integrate(lattice, A, B)
				Z3 = integrate(lattice, B, A)

				@test abs((Z1-Z2) / Z1) <= tol
				@test abs((Z1-Z3) / Z1) <= tol

				C = randomgmps(Float64, length(lattice), D=2)
				Z1 = integrate(lattice, A * B * C)
				Z2 = integrate(lattice, A, B, C)
				Z3 = integrate(lattice, B, A, C)

				@test abs((Z1-Z2) / Z1) <= tol
				@test abs((Z1-Z3) / Z1) <= tol
			end
		end		
	end

	trunc = truncdimcutoff(D=100, ϵ=1.0e-9)
	for N in (2,3)
		for bands in (1,2,3)
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt=0.2, bands=bands, contour=:real, ordering=ordering)
				A = randomgmps(ComplexF64, length(lattice), D=4)
				B = randomgmps(ComplexF64, length(lattice), D=6)
				C = mult(A, B, trunc=trunc)
				Z1 = integrate(lattice, C)
				Z2 = integrate(lattice, A, B)
				Z3 = integrate(lattice, B, A)

				@test abs((Z1-Z2) / Z1) <= tol
				@test abs((Z1-Z3) / Z1) <= tol

				C = randomgmps(ComplexF64, length(lattice), D=2)

				Z1 = integrate(lattice, A * B * C)
				Z2 = integrate(lattice, A, B, C)
				Z3 = integrate(lattice, B, C, A)

				@test abs((Z1-Z2) / Z1) <= tol
				@test abs((Z1-Z3) / Z1) <= tol
			end
		end		
	end	
end

@testset "GrassmannMPS: integration 4" begin
	rtol = 1.0e-4

	for N in (1, 2,3)
		for bands in (1,2)
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=2)
				B = randomgmps(Float64, length(lattice), D=4)
				C = randomgmps(Float64, length(lattice), D=4)
				D = randomgmps(Float64, length(lattice), D=2)

				AB = A * B 
				CD = C * D 

				Z1 = integrate(lattice, AB, CD)
				Z2 = integrate(lattice, A, B, C, D)

				@test abs((Z1-Z2) / Z1) <= rtol
			end
		end		
	end

	for N in (1,2)
		for bands in (1,2)
			for ordering in real_grassmann_orderings
				lattice = GrassmannLattice(N=N, δt=0.2, bands=bands, contour=:real, ordering=ordering)
				A = randomgmps(ComplexF64, length(lattice), D=2)
				B = randomgmps(ComplexF64, length(lattice), D=2)
				C = randomgmps(ComplexF64, length(lattice), D=2)
				D = randomgmps(ComplexF64, length(lattice), D=2)

				AB = A * B 
				CD = C * D 

				Z1 = integrate(lattice, AB, CD)
				Z2 = integrate(lattice, A, B, C, D)

				@test abs((Z1-Z2) / Z1) <= rtol
			end
		end		
	end	
end

# @testset "GrassmannMPS: integration 6" begin
# 	rtol = 1.0e-4

# 	for N in (1, 2,3)
# 		for bands in (1,2)
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
# 				Z2 = integrate(lattice, C, B, A, D, F, E)

# 				@test abs((Z1-Z2) / Z1) <= rtol
# 			end
# 		end		
# 	end

# 	for N in (1,2)
# 		for bands in (1,2)
# 			for ordering in real_grassmann_orderings
# 				lattice = GrassmannLattice(N=N, δt=0.2, bands=bands, contour=:real, ordering=ordering)
# 				A = randomgmps(ComplexF64, length(lattice), D=2)
# 				B = randomgmps(ComplexF64, length(lattice), D=2)
# 				C = randomgmps(ComplexF64, length(lattice), D=2)
# 				D = randomgmps(ComplexF64, length(lattice), D=2)
# 				E = randomgmps(ComplexF64, length(lattice), D=2)
# 				F = randomgmps(ComplexF64, length(lattice), D=2)


# 				ABC = A * B * C
# 				DEF = D * E * F 

# 				Z1 = integrate(lattice, ABC, DEF)
# 				Z2 = integrate(lattice, A, E, F, B, D, C)

# 				@test abs((Z1-Z2) / Z1) <= rtol
# 			end
# 		end		
# 	end	
# end


# @testset "GrassmannMPS: integration 7" begin
# 	rtol = 1.0e-4

# 	for N in (1, 2,3)
# 		for bands in (1,2)
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
# 				Z3 = integrate(lattice, D, B, C, A, G, F, E)

# 				@test abs((Z1-Z2) / Z1) <= rtol
# 				@test abs((Z1-Z3) / Z1) <= rtol
# 			end
# 		end		
# 	end

# 	for N in (1,2)
# 		for bands in (1,2)
# 			for ordering in real_grassmann_orderings
# 				lattice = GrassmannLattice(N=N, δt=0.2, bands=bands, contour=:real, ordering=ordering)
# 				A = randomgmps(ComplexF64, length(lattice), D=2)
# 				B = randomgmps(ComplexF64, length(lattice), D=2)
# 				C = randomgmps(ComplexF64, length(lattice), D=2)
# 				D = randomgmps(ComplexF64, length(lattice), D=2)
# 				E = randomgmps(ComplexF64, length(lattice), D=2)
# 				F = randomgmps(ComplexF64, length(lattice), D=2)
# 				G = randomgmps(ComplexF64, length(lattice), D=2)


# 				ABC = A * B * C
# 				DEFG = D * E * F * G

# 				Z1 = integrate(lattice, ABC, DEFG)
# 				Z2 = integrate(lattice, A, B, C, D, E, F, G)
# 				Z3 = integrate(lattice, D, B, C, A, G, F, E)

# 				@test abs((Z1-Z2) / Z1) <= rtol
# 				@test abs((Z1-Z3) / Z1) <= rtol
# 			end
# 		end		
# 	end	
# end
