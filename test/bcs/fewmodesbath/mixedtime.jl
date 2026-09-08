@testset "BCS, mixed time" begin
	U = 1.0; ϵ_d = 0.67; ω = 1.0; α = 0.5
	δt = 0.05; Nt = 5
	δτ = 0.1; Nτ = 5
	β = Nτ * δτ
	t = Nt * δt
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
	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
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
	@test relerr(g1, g1′) < 1.0e-2
	@test relerr(g2, g2′) < 5.0e-2

	# real and complex gaps: fast vs naive, then ED comparison (full thermal ρ)
	for Δ in (0.7, 0.3 + 0.4im)
		mpsK = (Δ isa Complex) ? complex(mpsK) : mpsK
		bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
		corr = correlationfunction(bath2, lattice)
		mpsI = hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc2)
		mpsI′ = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
		@test distance(mpsI, mpsI′) / norm(mpsI) < 1.0e-5
		cache = environments(lattice, mpsK, mpsI)
		g1 = [-im * cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.kt]
		g2 = [im * cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.kt]

		H, a, adag, H0 = bcs_ed(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)
		gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
		@test relerr(g1, gt_ed) < 3.0e-2
		@test relerr(g2, lt_ed) < 3.0e-2
	end
end
