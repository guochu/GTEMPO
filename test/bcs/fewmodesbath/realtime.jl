@testset "BCS, real time" begin
	U = 1.0; ϵ_d = 0.8; ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 5
	ts = [i*δt for i in 0:Nt]
	chi = 50; chi2 = 120
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10)

	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:Keldysh, bands=2)
	model = AndersonIM(U=U, μ=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:2
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)

	# trivial case: Δ = 0 reproduces two independent normal bands
	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
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

	# real and complex gaps: fast vs naive, then ED comparison
	for Δ in (0.6, 0.3 + 0.4im)
		mpsK = (Δ isa Complex) ? complex(mpsK) : mpsK
		bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
		corr = correlationfunction(bath2, lattice)
		mpsI = hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc2)
		mpsI′ = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
		@test distance(mpsI, mpsI′) / norm(mpsI) < 1.0e-5
		cache = environments(lattice, mpsK, mpsI)
		g1 = [-im * cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
		g2 = [-im * cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]

		H, a, adag, H0 = bcs_ed(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)
		gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)
		@test relerr(g1, gt_ed) < 3.0e-2
		@test relerr(g2, lt_ed) < 3.0e-2
	end
end
