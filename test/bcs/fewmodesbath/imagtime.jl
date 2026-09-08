@testset "BCS, imaginary time" begin
	U = 1.0; ϵ_d = 0.8; ω = 1.0; α = 0.5
	δτ = 0.1; Nτ = 10; β = Nτ * δτ
	chi = 60; chi2 = 200
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10)

	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
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
	g1 = gτ_series(lattice, mpsK, mpsI1, mpsI2)

	bath0 = bcsbath(bath, Δ=0)
	corr = correlationfunction(bath0, lattice)
	mpsI = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc2)
	g2 = gτ_series(lattice, mpsK, mpsI)
	@test relerr(g1, g2) < 1.0e-2

	# real gap: naive and fast IF constructions agree, compare with ED
	for (Δ, tag) in ((0.6, "real gap"), (0.2 + 0.45im, "complex gap"))
		mpsK = (Δ isa Complex) ? complex(mpsK) : mpsK
		bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
		corr = correlationfunction(bath2, lattice)
		mpsI = hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc2)
		mpsI′ = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
		@test distance(mpsI, mpsI′) / norm(mpsI) < 1.0e-5
		cache = environments(lattice, mpsK, mpsI)
		g1 = cached_gf_fast(lattice, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)

		H, a, adag, H0 = bcs_ed(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)
		g2 = gτ_ed(H, a, adag, 0:δτ:β, β)
		@test relerr(g1, g2) < 1.0e-2
	end

	# ordering dependence: BCS IF on the alternative ordering
	lattice2 = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2, ordering=imag_orderings[2])
	Δ = 0.6
	bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
	corr2 = correlationfunction(bath2, lattice2)
	mpsK2 = sysdynamics(lattice2, model, trunc=trunc)
	for band in 1:2
		mpsK2 = boundarycondition!(mpsK2, lattice2, band=band, trunc=trunc)
	end
	mpsI2 = hybriddynamics(lattice2, corr2, orbital=1, trunc=trunc)
	g3 = gτ_series(lattice2, mpsK2, mpsI2)
	H, a, adag, H0 = bcs_ed(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)
	@test relerr(g3, gτ_ed(H, a, adag, 0:δτ:β, β)) < 1.0e-2
end
