@testset "BCS (continuous bath), real time" begin
	U = 1.0; ϵ_d = 0.8
	β = 1.0; δt = 0.05; Nt = 5
	ts = [i*δt for i in 0:Nt]
	chi = 50; chi2 = 120
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10)

	# trivial case: Δ = 0 reproduces two independent normal bands
	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:Keldysh, bands=2)
	model = AndersonIM(U=U, μ=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:2
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)
	bath = fermionicbath(semicircular(t=1), β=β)
	corr = correlationfunction(bath, lattice)
	mpsI1 = hybriddynamics(lattice, corr, band=1, trunc=trunc)
	mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)
	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	g1 = [-im * cached_greater(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.k]
	g2 = [-im * cached_lesser(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.k]

	bath0 = bcsbath(bath, Δ=0)
	corr = correlationfunction(bath0, lattice)
	mpsI = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc2)
	cache = environments(lattice, mpsK, mpsI)
	g1′ = [-im * cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	g2′ = [-im * cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	@test relerr(g1, g1′) < 1.0e-2
	@test relerr(g2, g2′) < 5.0e-2

	# quadratic model with a complex gap: Toulouse ED on the discretized bath
	# (empty initial state, nsys = 0)
	Δ = 0.2 + 0.45im
	model = AndersonIM(U=0, μ=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:2
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	bath2 = bcsbath(fermionicbath(semicircular(t=1), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics(lattice, corr, orbital=1, trunc=trunc2)
	cache = environments(lattice, mpsK, mpsI)
	g1 = [-im * cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	g2 = [-im * cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]

	disbath2 = discretebath(bath2, δw=0.02)
	ed_model = Toulouse(disbath2, ϵ_d=-ϵ_d)
	g1′, g2′ = toulouse_neq_greater_lesser(ed_model, ts, nsys=0)
	@test relerr(g1, g1′) < 2.0e-2
	@test norm(g2 - g2′) < 5.0e-2
end
