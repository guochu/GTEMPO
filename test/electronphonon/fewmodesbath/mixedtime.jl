@testset "Electron-phonon, mixed time: independent bosons (phonon only)" begin
	μ = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=120, ϵ=1.0e-10)
	spec = DiracDelta(ω=1, α=0.5)

	for (U, bands) in ((0.0, 1), (1.0, 2))
		lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=bands)
		flat = FockLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=bands)
		bath = bosonicbath(spec, β=β)
		corr = correlationfunction(bath, flat)
		mpsI = hybriddynamics(flat, corr, trunc=trunc)

		model = (U == 0) ? ToulouseIM(μ=μ) : AndersonIM(U=U, μ=μ)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		adt = reweighting!(lattice, mpsK, flat, mpsI, trunc=trunc)
		for band in 1:bands
			adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
		end

		cache = environments(lattice, adt)
		gt = [-im * cached_gf(lattice, (ContourIndex(k, conj=false, branch=:+, band=1), ContourIndex(1, conj=true, branch=:+, band=1)), adt; cache=cache) for k in 1:Nt+1]
		lt = [im * cached_gf(lattice, (ContourIndex(1, conj=true, branch=:-, band=1), ContourIndex(k, conj=false, branch=:+, band=1)), adt; cache=cache) for k in 1:Nt+1]
		gτ = [cached_gf(lattice, (ContourIndex(k, conj=false, branch=:τ, band=1), ContourIndex(1, conj=true, branch=:τ, band=1)), adt; cache=cache) for k in 1:Nτ+1]

		gt_ed = [independentbosons_greater_μ(spec, tj, β=β, μ=μ, U=U, bands=bands) for tj in ts]
		lt_ed = [independentbosons_lesser_μ(spec, tj, β=β, μ=μ, U=U, bands=bands) for tj in ts]
		gτ_ed = independentbosons_Gτ_μ(spec, β=β, μ=μ, Nτ=Nτ, U=U, bands=bands)

		@test relerr(gt, gt_ed) < 3.0e-2
		@test relerr(lt, lt_ed) < 3.0e-2
		@test relerr(gτ, gτ_ed) < 1.0e-2
	end
end

@testset "Electron-phonon, mixed time: fermionic + phonon bath" begin
	μ = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=120, ϵ=1.0e-10)
	ωf, αf = 1.0, 0.5
	ω₀, α = 1.0, 0.5

	for (U, bands) in ((0.0, 1), (1.0, 2))
		# ED reference with the interacting thermal state
		H, a, adag, H0 = mixedbath_ed(μ=μ, U=U, ωf=ωf, αf=αf, ω₀=ω₀, α=α, d=8)
		gt_ed, lt_ed = greater_lesser_ed_mixed(H, a, adag, ts, β)
		gτ_ref = gτ_ed(H, a, adag, 0:δτ:β, β)

		lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=bands)
		flat = FockLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=bands)

		fbath = fermionicbath(DiracDelta(ω=ωf, α=αf), β=β)
		lattice1 = similar(lattice, bands=1)
		fcorr = correlationfunction(fbath, lattice1)
		mpsI_e = hybriddynamics(lattice1, fcorr, trunc=trunc)
		Is = [fillband(lattice, mpsI_e, band=b) for b in 1:bands]

		pbath = bosonicbath(DiracDelta(ω=ω₀, α=α), β=β)
		pcorr = correlationfunction(pbath, flat)
		mpsI_p = hybriddynamics(flat, pcorr, trunc=trunc)

		model = (U == 0) ? ToulouseIM(μ=μ) : AndersonIM(U=U, μ=μ)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		adt = reweighting!(lattice, mpsK, flat, mpsI_p, trunc=trunc)
		for band in 1:bands
			adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
		end

		cache = environments(lattice, adt, Is...)
		gt = [-im * cached_gf(lattice, (ContourIndex(k, conj=false, branch=:+, band=1), ContourIndex(1, conj=true, branch=:+, band=1)), adt, Is...; cache=cache) for k in 1:Nt+1]
		lt = [im * cached_gf(lattice, (ContourIndex(1, conj=true, branch=:-, band=1), ContourIndex(k, conj=false, branch=:+, band=1)), adt, Is...; cache=cache) for k in 1:Nt+1]
		gτ = cached_gf_fast(lattice, adt, Is...; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
		gτ[end] = 1 - gτ[1]

		@test relerr(gt, gt_ed) < 3.0e-2
		@test relerr(lt, lt_ed) < 3.0e-2
		# d = 8 phonon truncation in the ED reference limits the accuracy
		@test relerr(gτ, gτ_ref) < 2.0e-2
	end
end
