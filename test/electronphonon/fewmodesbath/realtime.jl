@testset "Electron-phonon, real time: independent bosons (phonon only)" begin
	μ = 0.5
	β = 1.0; δt = 0.05; Nt = 10
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)
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

		model = AndersonIM(U=U, μ=μ)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)
		for band in 1:bands
			mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
		end
		adt = reweighting!(lattice, mpsK, flat, mpsI, trunc=trunc)

		cache = environments(lattice, adt)
		gt = [-im * cached_greater(lattice, k, adt, band=1, cache=cache) for k in 1:Nt+1]
		lt = [-im * cached_lesser(lattice, k, adt, band=1, cache=cache) for k in 1:Nt+1]
		@test relerr(gt, gt_ed) < 3.0e-2
		@test relerr(lt, lt_ed) < 3.0e-2
	end
end

@testset "Electron-phonon, real time: fermionic + phonon bath" begin
	μ = 0.5
	β = 1.0; δt = 0.05; Nt = 10
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)
	ωf, αf = 1.0, 0.5
	ω₀, α = 1.0, 0.5

	for (U, bands) in ((0.0, 1), (1.0, 2))
		# ED reference for the combined model
		H, a, adag, H0 = mixedbath_ed(μ=μ, U=U, ωf=ωf, αf=αf, ω₀=ω₀, α=α, d=8)
		gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

		lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1, bands=bands)
		flat = FockLattice(N=Nt, δt=δt, contour=:real, order=1, bands=bands)

		fbath = fermionicbath(DiracDelta(ω=ωf, α=αf), β=β)
		lattice1 = similar(lattice, bands=1)
		fcorr = correlationfunction(fbath, lattice1)
		mpsI_e = hybriddynamics(lattice1, fcorr, trunc=trunc)
		Is = [fillband(lattice, mpsI_e, band=b) for b in 1:bands]

		pbath = bosonicbath(DiracDelta(ω=ω₀, α=α), β=β)
		pcorr = correlationfunction(pbath, flat)
		mpsI_p = hybriddynamics(flat, pcorr, trunc=trunc)

		model = AndersonIM(U=U, μ=μ)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)
		for band in 1:bands
			mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
		end
		adt = reweighting!(lattice, mpsK, flat, mpsI_p, trunc=trunc)

		cache = environments(lattice, adt, Is...)
		gt = [-im * cached_greater(lattice, k, adt, Is..., band=1, cache=cache) for k in 1:Nt+1]
		lt = [-im * cached_lesser(lattice, k, adt, Is..., band=1, cache=cache) for k in 1:Nt+1]
		@test relerr(gt, gt_ed) < 3.0e-2
		@test relerr(lt, lt_ed) < 3.0e-2
	end
end
