println("------------------------------------")
println("|            fillband              |")
println("------------------------------------")



@testset "IF-imaginary time: fillband" begin
	N = 25
	δτ = 0.01
	ϵ_d = 1.25*pi
	β = N * δτ

	rtol = 1.0e-4

	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)

	bath = fermionicbath(spectrum_func(), β=β, μ=1)
	corr = Δτ(bath, N=N, δτ=δτ)
	for ordering in imag_grassmann_orderings
		lattice = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=ordering, bands=3)
		lattice2 = similar(lattice, bands=1)
		mpsI2 = hybriddynamics(lattice2, corr, trunc=trunc) 
		for band in 1:lattice.bands
			mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=band) 
			mpsI1′ = fillband(lattice, mpsI2, band=band)
			@test distance(mpsI1, mpsI1′) / norm(mpsI1) < rtol
		end
	end
end


@testset "IF-real time: fillband" begin
	N = 10
	δt = 0.01
	t = N * δt
	ϵ_d = 1.25*pi
	β = 1.
	rtol = 1.0e-4

	trunc = truncdimcutoff(D=100, ϵ=1.0e-10, add_back=0)
		
	# println("μ = ", μ)
	bath = fermionicbath(spectrum_func(), β=β, μ=0.)
	corr = Δt(bath, N=N, t=t)

	for ordering in real_ac_grassmann_orderings
		lattice = GrassmannLattice(N=N, δt=δt, contour=:real, ordering=ordering, bands=2)
		lattice2 = similar(lattice, bands=1)
		mpsI2 = hybriddynamics(lattice2, corr, trunc=trunc) 
		for band in 1:lattice.bands
			mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=band) 
			mpsI1′ = fillband(lattice, mpsI2, band=band)
			@test distance(mpsI1, mpsI1′) / norm(mpsI1) < rtol
		end
	end
end


@testset "IF-mixed time: fillband" begin

	N = 25
	δτ = 0.01
	ϵ_d = 1.25*pi
	dw = 0.1
	β = N * δτ
	δt = 0.01
	Nt = 20
	t = δt * Nt

	rtol = 1.0e-4

	trunc = truncdimcutoff(D=100, ϵ=1.0e-10, add_back=0)
		
	bath = fermionicbath(spectrum_func(), β=β, μ=0)

	exact_model = AndersonIM(μ=ϵ_d, U=1)
	corr = Δm(bath, Nτ=N, Nt=Nt, t=t)
	for ordering in mixed_ac_grassmann_orderings
		lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=N, δτ=β/N, contour=:mixed, ordering=ordering, bands=2)
		lattice2 = similar(lattice, bands=1)
		mpsI2 = hybriddynamics(lattice2, corr, trunc=trunc) 
		for band in 1:lattice.bands
			mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=band) 
			mpsI1′ = fillband(lattice, mpsI2, band=band)
			@test distance(mpsI1, mpsI1′) / norm(mpsI1) < rtol
		end		
	end
end