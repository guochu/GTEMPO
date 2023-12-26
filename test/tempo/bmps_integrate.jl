println("------------------------------------")
println("|          BMPS integrate          |")
println("------------------------------------")



@testset "GrassmannMPS: BMPS integration" begin
	rtol = 1.0e-4

	for N in (1, 2,3)
		for bands in (1,2)
			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=0.1, bands=bands, contour=:imag, ordering=ordering)
				A = randomgmps(Float64, length(lattice), D=2)
				B = randomgmps(Float64, length(lattice), D=4)
				C = randomgmps(Float64, length(lattice), D=4)
				D = randomgmps(Float64, length(lattice), D=2)

				Z1 = integrate(lattice, A, B)
				Z2 = integrate(lattice, A, B, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol

				Z1 = integrate(lattice, A, B, C)
				Z2 = integrate(lattice, A, B, C, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol		

				Z1 = integrate(lattice, A, B, C, D)
				Z2 = integrate(lattice, A, B, C, D, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol		

				Z1 = integrate(lattice, [A, B], D)
				Z2 = integrate(lattice, [A, B], D, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol			

				Z1 = integrate(lattice, [A, B, C], D)
				Z2 = integrate(lattice, [A, B, C], D, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol	

				Z1 = integrate(lattice, [A, B, C, D], D, A)
				Z2 = integrate(lattice, [A, B, C, D], D, A, alg=BMPSIntegrate())
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

				Z1 = integrate(lattice, A, B)
				Z2 = integrate(lattice, A, B, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol

				Z1 = integrate(lattice, A, B, C)
				Z2 = integrate(lattice, A, B, C, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol		

				Z1 = integrate(lattice, A, B, C, D)
				Z2 = integrate(lattice, A, B, C, D, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol		

				Z1 = integrate(lattice, [A, B], D)
				Z2 = integrate(lattice, [A, B], D, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol			

				Z1 = integrate(lattice, [A, B, C], D)
				Z2 = integrate(lattice, [A, B, C], D, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol	

				Z1 = integrate(lattice, [A, B, C, D], D, A)
				Z2 = integrate(lattice, [A, B, C, D], D, A, alg=BMPSIntegrate())
				@test abs((Z1-Z2) / Z1) <= rtol																									
			end
		end		
	end	
end