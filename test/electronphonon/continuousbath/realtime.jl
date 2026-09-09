@testset "Electron-phonon (continuous bath), real time: independent bosons" begin
	μ = 0.5
	β = 1.0; δt = 0.05; Nt = 10
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)
	# The Keldysh contour starts from a separable thermal state, while the
	# analytic equilibrium solution of a continuous spectrum assumes the
	# interacting equilibrium — for the real-time check we therefore use a
	# δ-mode phonon with an exact ED reference (the continuous Leggett
	# spectrum is covered by the imaginary- and mixed-time tests below).
	spec = DiracDelta(ω=1, α=0.5)

	for (U, bands) in ((0.0, 1), (1.0, 2))
		# ED reference (phonon truncated to d = 8 levels)
		H, a, adag, H0 = phonon_ed(μ=μ, U=U, ω₀=1.0, α=0.5, d=8)
		gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

		lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1, bands=bands)
		flat = FockLattice(N=Nt, δt=δt, contour=:real, order=1, bands=bands)
		bath = bosonicbath(spec, β=β)
		corr = correlationfunction(bath, flat)
		mpsI = hybriddynamics(flat, corr, trunc=trunc)

		model = (U == 0) ? ToulouseIM(μ=μ) : AndersonIM(U=U, μ=μ)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)
		for band in 1:bands
			mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
		end
		adt = reweighting!(lattice, mpsK, flat, mpsI, trunc=trunc)

		cache = environments(lattice, adt)
		gt = [-im * cached_greater(lattice, k, adt, band=1, cache=cache) for k in 1:Nt+1]
		lt = [-im * cached_lesser(lattice, k, adt, band=1, cache=cache) for k in 1:Nt+1]
		@test relerr(gt, gt_ed) < 5.0e-2
		@test relerr(lt, lt_ed) < 5.0e-2
	end
end
