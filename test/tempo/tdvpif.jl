println("------------------------------------")
println("|              TDVPIF              |")
println("------------------------------------")

# TDVPIF views the influence functional as the "equilibrium state" IF = exp(H)
# of the influence operator H and computes it with a second-order single-site
# TDVP imaginary-time flow dz/dτ = H·z from τ = 0 (identity IF) to τ = 1.

@testset "TDVPIF: imaginary-time" begin

	N = 5
	δτ = 0.1
	β = N * δτ

	rtol = 5.0e-2

	trunc = truncdimcutoff(D=100, ϵ=1.0e-6, add_back=0)

	base_alg = PartialIF(trunc=trunc)

	for μ in (0, 5)
		for spec in (spectrum_func(), spectrum_func2())

			bath = fermionicbath(spec, β=β, μ=μ)

			for ordering in imag_grassmann_orderings
				lattice = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=ordering)

				corr = correlationfunction(bath, lattice)
				mpsI1 = hybriddynamics(lattice, corr, base_alg)

				for δ in (0.1, 0.05)
					mpsI2 = hybriddynamics(lattice, corr, TDVPIF(trunc=trunc, δ=δ))
					@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
				end
			end
		end
	end
end

@testset "TDVPIF: real-time" begin

	N = 3
	δt = 0.1
	β = 1

	rtol = 5.0e-2
	trunc = truncdimcutoff(D=48, ϵ=1.0e-8)

	base_alg = PartialIF(trunc=trunc)

	# the real-time influence operator is only defined for these orderings
	tdvpif_orderings = [A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1a1ā1b1b̄1()]

	for spec in (spectrum_func(1), spectrum_func2(1))

		bath = fermionicbath(spec, β=β, μ=0)

		for ordering in tdvpif_orderings
			lattice = GrassmannLattice(N=N, δt=δt, contour=:real, ordering=ordering)

			corr = correlationfunction(bath, lattice)
			mpsI1 = hybriddynamics(lattice, corr, base_alg)

			for δ in (0.1, 0.05)
				mpsI2 = hybriddynamics(lattice, corr, TDVPIF(trunc=trunc, δ=δ))
				@test distance(mpsI1, mpsI2) / norm(mpsI1) < rtol
			end
		end
	end
end

@testset "TDVPIF: convergence with δ" begin

	# the flow error decreases with the step size δ; the residual against
	# ExactTTIIF (exact exponentiation of the same influence operator) is
	# the pure TDVP integrator error
	N = 5
	β = 0.5
	trunc = truncdimcutoff(D=64, ϵ=1.0e-8)

	bath = fermionicbath(spectrum_func(), β=β, μ=0)
	lattice = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=A1Ā1B1B̄1())
	corr = correlationfunction(bath, lattice)

	mpsE = hybriddynamics(lattice, corr, ExactTTIIF(algmult=SVDCompression(trunc)))
	errs = [distance(hybriddynamics(lattice, corr, TDVPIF(trunc=trunc, δ=δ)), mpsE) / norm(mpsE) for δ in (0.2, 0.1, 0.05)]
	@test errs[1] > errs[2] > errs[3]
	@test errs[3] < 1.0e-3
end

@testset "TDVPIF: in-place flow" begin

	# the flow evolves z(τ=1) = e^H·z(0) in place: merging the influence
	# functional into an arbitrary initial state in a single flow must agree
	# with multiplying the separately constructed IF onto it
	N = 3
	β = 1
	trunc = truncdimcutoff(D=48, ϵ=1.0e-8)

	bath = fermionicbath(spectrum_func(1), β=β, μ=0)
	lattice = GrassmannLattice(N=N, δt=0.1, contour=:real, ordering=A1Ā1a1ā1B1B̄1b1b̄1())
	corr = correlationfunction(bath, lattice)

	Random.seed!(12354)
	z0 = randomgmps(ComplexF64, length(lattice), D=2)
	mpsI = hybriddynamics(lattice, corr, TDVPIF(trunc=trunc, δ=0.05))

	z1 = hybriddynamics!(copy(z0), lattice, corr, TDVPIF(trunc=trunc, δ=0.05))
	z2 = mult(mpsI, z0, SVDCompression(trunc))
	@test distance(z1, z2) / norm(z2) < 1.0e-2
end
