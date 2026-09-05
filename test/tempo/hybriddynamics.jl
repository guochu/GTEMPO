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
	algs = [XTRGIF(k=5, fast=true), XTRGIF(k=5, fast=false), ExactTTIIF()]
		
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
	alg2 = XTRGIF(k=5, algevo=WII(), algmult=SVDCompression(trunc))
	alg3 = XTRGIF(k=5, algmult=DMRG1(trunc=trunc, initguess=:svd))
	alg4 = XTRGIF(k=5, algmult=DMRG1(trunc=trunc, initguess=:pre))
	alg5 = XTRGIF(k=5, algmult=DMRG1(trunc=trunc, initguess=:rand, maxiter=10))
	alg6 = XTRGIF(k=5, algmult=DMRG1(trunc=trunc), fast=false)
	alg7 = XTRGIF(k=5, algevo=ComplexStepper(WII()), algmult=DMRG2(trunc=trunc, initguess=:svd))
	alg8 = ExactTTIIF(algmult=DMRG1(trunc=trunc, initguess=:rand))
	alg9 = ExactTTIIF(algmult=SVDCompression(trunc))

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