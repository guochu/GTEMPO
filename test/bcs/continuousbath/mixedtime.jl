@testset "BCS (continuous bath), mixed time" begin
	U = 1.0; ϵ_d = 0.67
	δt = 0.05; Nt = 5
	δτ = 0.1; Nτ = 5
	β = Nτ * δτ
	ts = [i*δt for i in 0:Nt]
	chi = 50; chi2 = 120
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10)

	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:Kadanoff, bands=2)
	model = AndersonIM(U=U, μ=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:2
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end

	# trivial case: Δ = 0 reproduces two independent normal bands
	bath = fermionicbath(semicircular(t=1), β=β)
	corr = correlationfunction(bath, lattice)
	mpsI1 = hybriddynamics(lattice, corr, band=1, trunc=trunc)
	mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)
	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	g1 = [-im * cached_greater(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.kt]
	g2 = [im * cached_lesser(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.kt]

	bath0 = bcsbath(bath, Δ=0)
	corr = correlationfunction(bath0, lattice)
	mpsI = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc2)
	cache = environments(lattice, mpsK, mpsI)
	g1′ = [-im * cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.kt]
	g2′ = [im * cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.kt]
	@test relerr(g1, g1′) < 2.0e-2
	@test relerr(g2, g2′) < 5.0e-2

	# quadratic model with a complex gap: Toulouse ED on the discretized bath
	Δ = 0.1 + 0.2im
	model = AndersonIM(U=0, μ=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:2
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	bath2 = bcsbath(fermionicbath(semicircular(t=1), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
	cache = environments(lattice, mpsK, mpsI)
	g1 = [-im * cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.kt]
	g2 = [im * cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.kt]

	disbath2 = discretebath(bath2, δw=0.02)
	ed_model = Toulouse(disbath2, ϵ_d=-ϵ_d)
	g1_ref, g2_ref = toulouse_greater_lesser(ed_model, ts)
	g1_ref, g2_ref = im * g1_ref, -im * g2_ref
	# convention: the equilibrium toulouse greater/lesser equal the raw
	# cached_greater/cached_lesser correlators times (im, -im)
	g1 = [cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.kt]
	g2 = [cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.kt]
	# the δw discretization of the bath dominates the residual error
	@test relerr(g1, g1_ref) < 4.0e-2
	@test relerr(g2, g2_ref) < 2.0e-2
end
