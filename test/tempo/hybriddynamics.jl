println("------------------------------------")
println("|         Hybrid dynamics          |")
println("------------------------------------")


@testset "Hybriddynamics: imaginary-time" begin

	N = 6
	δτ = 0.1
	β = N * δτ
	rtol = 1.0e-2

	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

	base_alg = PartialIF(trunc=trunc)
	algs = [TranslationInvariantIF(k=5, fast=true), TranslationInvariantIF(k=5, fast=false), ExactTranslationInvariantIF()]
		
	for μ in (-5, 0, 5)
		# println("μ = ", μ)
		for spec in (spectrum_func(), spectrum_func2())

			bath = fermionicbath(spec, β=β, μ=μ)

			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=ordering)

				corr = correlationfunction(bath, lattice)
				mpsI1 = hybriddynamics(lattice, corr, base_alg) 

				mpsI0 = hybriddynamics_naive(lattice, corr, trunc=trunc) 
				@test distance(mpsI1, mpsI0) / norm(mpsI1) < rtol


				for alg in algs
					mpsI2 = hybriddynamics(lattice, corr, alg)
					@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
				end
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

	base_alg = PartialIF(trunc=trunc)
	alg2 = TranslationInvariantIF(k=5, algevo=WII(), algmult=SVDCompression(trunc))
	alg3 = TranslationInvariantIF(k=5, algmult=DMRGMult1(trunc=trunc, initguess=:svd))
	alg4 = TranslationInvariantIF(k=5, algmult=DMRGMult1(trunc=trunc, initguess=:pre))
	alg5 = TranslationInvariantIF(k=5, algmult=DMRGMult1(trunc=trunc, initguess=:rand, maxiter=10))
	alg6 = TranslationInvariantIF(k=5, algmult=DMRGMult1(trunc=trunc), fast=false)
	alg7 = TranslationInvariantIF(k=5, algevo=ComplexStepper(WII()), algmult=DMRGMult2(trunc=trunc, initguess=:svd))
	alg8 = ExactTranslationInvariantIF(algmult=DMRGMult1(trunc=trunc, initguess=:rand))
	alg9 = ExactTranslationInvariantIF(algmult=SVDCompression(trunc))

	algs = [alg2, alg3, alg4, alg5, alg6, alg7, alg8, alg9]

		

	for spec in (spectrum_func(1), spectrum_func2(1))

		bath = fermionicbath(spec, β=β, μ=0)

		for ordering in real_grassmann_orderings
			# println("ordering is ", ordering)
			lattice = GrassmannLattice(N=N, δt=δt, contour=:real, ordering=ordering)

			corr = correlationfunction(bath, lattice)
			mpsI1 = hybriddynamics(lattice, corr, base_alg) 

			mpsI0 = hybriddynamics_naive(lattice, corr, trunc=trunc) 
			@test distance(mpsI1, mpsI0) / norm(mpsI1) < rtol

			for alg in algs
				mpsI2 = hybriddynamics(lattice, corr, alg)
				@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
			end

		end
	end

end