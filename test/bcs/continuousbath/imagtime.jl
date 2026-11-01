@testset "BCS (continuous bath), imaginary time" begin
	U = 1.0; ϵ_d = 0.8
	δτ = 0.1; Nτ = 10; β = Nτ * δτ
	chi = 60; chi2 = 200
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10)

	# trivial case: Δ = 0 reproduces two independent normal bands
	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
	model = AndersonIM(U=U, ϵ_d=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:2
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	bath = fermionicbath(semicircular(t=1), β=β)
	corr = correlationfunction(bath, lattice)
	mpsI1 = hybriddynamics(lattice, corr, band=1, trunc=trunc)
	mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)
	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	g1 = cached_gf_fast(lattice, mpsK, mpsI1, mpsI2; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)

	bath0 = bcsbath(bath, Δ=0)
	corr = correlationfunction(bath0, lattice)
	mpsI = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc2)
	cache = environments(lattice, mpsK, mpsI)
	g2 = cached_gf_fast(lattice, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
	# same BCS τ vertex-shift residual as the few-mode case: the two channels
	# agree up to the (slightly different) O(δτ) discretization constants
	@test relerr(g1, g2) < 2.0e-2

	# quadratic model with a finite complex gap: compare with the Toulouse ED
	# on the discretized BCS bath, for both orderings
	Δ = 0.2 + 0.45im
	bath2 = bcsbath(fermionicbath(semicircular(t=1), β=β), Δ=Δ)
	gs = []
	for ordering in imag_orderings
		lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2, ordering=ordering)
		model = AndersonIM(U=0, ϵ_d=-ϵ_d)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		mpsK = complex(mpsK)
		for band in 1:2
			mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
		end
		corr = correlationfunction(bath2, lattice)
		mpsI = hybriddynamics(lattice, corr, orbital=1, trunc=trunc2)
		push!(gs, gτ_series(lattice, mpsK, mpsI))
	end
	disbath2 = discretebath(bath2, δw=0.02)
	ed_model = Toulouse(disbath2, ϵ_d=-ϵ_d)
	g_ref = toulouse_Gτ(ed_model, collect(0:δτ:β))
	# the former ≈2.5% (δw-insensitive) systematic residual was the τ-ring
	# VERTEX shift, not an anomalous-block convention issue: with the τ ring
	# reverted to the raw read (2026-10-30) it drops to 8e-4
	@test relerr(gs[1], g_ref) < 5.0e-3
	@test relerr(gs[2], g_ref) < 5.0e-3
end
