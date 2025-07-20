println("------------------------------------")
println("|       Mult integrateband         |")
println("------------------------------------")



@testset "IF-imaginary time: integrateband" begin
	N = 25
	δτ = 0.01
	ϵ_d = 1.25*pi
	β = N * δτ

	rtol = 1.0e-4

	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

	trunc2 = truncdimcutoff(D=300, ϵ=1.0e-10, add_back=0)
	algmult = DMRGMult1(trunc=trunc2, verbosity=1, maxiter=5)

	bath1 = fermionicbath(spectrum_func(), β=β, μ=0)
	corr1 = Δτ(bath1, N=N, δτ=δτ)
	bath2 = fermionicbath(spectrum_func(), β=β, μ=1)
	corr2 = Δτ(bath2, N=N, δτ=δτ)
	for ordering in imag_ac_grassmann_orderings
		lattice = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=ordering, bands=2)
		mpsI1 = hybriddynamics(lattice, corr1, trunc=trunc, band=1)
		mpsI2 = hybriddynamics(lattice, corr2, trunc=trunc, band=2)

		for band in 1:2
			res1 = integrateband(lattice, mult(mpsI1, mpsI2), band=band)
			res2 = integrateband(lattice, mpsI1, mpsI2, algmult; band=band)
			@test distance(res1, res2) / norm(res1) < rtol	
		end
	end
end

@testset "IF-real time: integrateband" begin
	N = 6
	δt = 0.01
	t = N * δt
	ϵ_d = 1.25*pi
	β = 1.
	rtol = 1.0e-4

	trunc = truncdimcutoff(D=100, ϵ=1.0e-10, add_back=0)
	
	trunc2 = truncdimcutoff(D=300, ϵ=1.0e-10, add_back=0)
	algmult = DMRGMult1(trunc=trunc2, verbosity=1, maxiter=5)

	bath1 = fermionicbath(spectrum_func(), β=β, μ=0.)
	corr1 = Δt(bath1, N=N, t=t)
	bath2 = fermionicbath(spectrum_func(), β=β, μ=1.)
	corr2 = Δt(bath2, N=N, t=t)

	for ordering in real_ac_grassmann_orderings
		lattice = GrassmannLattice(N=N, δt=δt, contour=:real, ordering=ordering, bands=2)
		mpsI1 = hybriddynamics(lattice, corr1, trunc=trunc, band=1)
		mpsI2 = hybriddynamics(lattice, corr2, trunc=trunc, band=2)

		for band in 1:2
			res1 = integrateband(lattice, mult(mpsI1, mpsI2), band=band)
			res2 = integrateband(lattice, mpsI1, mpsI2, algmult; band=band)
			@test distance(res1, res2) / norm(res1) < rtol	
		end
	end
end

@testset "IF-mixed time: integrateband" begin

	N = 6
	δτ = 0.01
	ϵ_d = 1.25*pi
	dw = 0.1
	β = N * δτ
	δt = 0.01
	Nt = 6
	t = δt * Nt

	rtol = 1.0e-4

	trunc = truncdimcutoff(D=100, ϵ=1.0e-10, add_back=0)
	
	trunc2 = truncdimcutoff(D=300, ϵ=1.0e-10, add_back=0)
	algmult = DMRGMult1(trunc=trunc2, verbosity=1, maxiter=5)

	bath1 = fermionicbath(spectrum_func(), β=β, μ=0)
	corr1 = Δm(bath1, Nτ=N, Nt=Nt, t=t)
	bath2 = fermionicbath(spectrum_func(), β=β, μ=1)
	corr2 = Δm(bath2, Nτ=N, Nt=Nt, t=t)

	for ordering in mixed_ac_grassmann_orderings
		lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=N, δτ=β/N, contour=:mixed, ordering=ordering, bands=2)
		mpsI1 = hybriddynamics(lattice, corr1, trunc=trunc, band=1)
		mpsI2 = hybriddynamics(lattice, corr2, trunc=trunc, band=2)

		for band in 1:2
			res1 = integrateband(lattice, mult(mpsI1, mpsI2), band=band)
			res2 = integrateband(lattice, mpsI1, mpsI2, algmult; band=band)
			@test distance(res1, res2) / norm(res1) < rtol	
		end
	end
end