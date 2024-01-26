println("------------------------------------")
println("|         Hybrid dynamics          |")
println("------------------------------------")


@testset "Hybriddynamics: imaginary-time" begin

	N = 6
	δτ = 0.1
	β = N * δτ
	rtol = 1.0e-2

	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

	alg1 = PartialIF(trunc=trunc)
	alg2 = TranslationInvariantIF(k=5)
		
	for μ in (-5, 0, 5)
		# println("μ = ", μ)
		for spec in (spectrum_func(), spectrum_func2())

			bath = fermionicbath(spec, β=β, μ=μ)

			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=ordering)

				corr = correlationfunction(bath, lattice)
				mpsI1 = hybriddynamics(lattice, corr, alg1) 
				mpsI2 = hybriddynamics(lattice, corr, alg2)

				@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
			end
		
		end
	end
end

@testset "Hybriddynamics: real-time" begin

	N = 3
	δt = 0.1
	β = 1

	rtol = 1.0e-2
	trunc = truncdimcutoff(D=100, ϵ=1.0e-6, add_back=0)

	alg1 = PartialIF(trunc=trunc)
	alg2 = TranslationInvariantIF(k=5, algevo=WII(), algmult=SVDCompression(trunc))
	alg3 = TranslationInvariantIF(k=5, algmult=DMRG1(D=trunc.D, tol=trunc.ϵ))
	alg4 = TranslationInvariantIF(k=5, algevo=ComplexStepper(WII()), algmult=DMRG2(trunc=trunc))

		

	for spec in (spectrum_func(1), spectrum_func2(1))

		bath = fermionicbath(spec, β=β, μ=0)

		for ordering in real_grassmann_orderings
			# println("ordering is ", ordering)
			lattice = GrassmannLattice(N=N, δt=δt, contour=:real, ordering=ordering)

			corr = correlationfunction(bath, lattice)
			mpsI1 = hybriddynamics(lattice, corr, alg1) 
			mpsI2 = hybriddynamics(lattice, corr, alg2)

			# println(norm(mpsI1), " ", norm(mpsI2), " ", distance(mpsI1, mpsI2))
			@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol

			mpsI3 = hybriddynamics(lattice, corr, alg3)
			@test distance(mpsI1, mpsI3) / norm(mpsI1) < rtol

			mpsI4 = hybriddynamics(lattice, corr, alg4)
			@test distance(mpsI1, mpsI4) / norm(mpsI1) < rtol
		end
	end

end