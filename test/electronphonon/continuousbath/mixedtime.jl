@testset "Electron-phonon (continuous bath), mixed time: independent bosons" begin
	μ = 0.5
	β = 1.0; δτ = 0.1; Nτ = round(Int, β/δτ)
	δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=120, ϵ=1.0e-10)
	spec = Leggett(d=3, ωc=1)

	for (U, bands) in ((0.0, 1), (1.0, 2))
		lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, bands=bands)
		flat = FockLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=bands)
		bath = bosonicbath(spec, β=β)
		corr = correlationfunction(bath, flat)
		mpsI = hybriddynamics(flat, corr, trunc=trunc)

		model = (U == 0) ? ToulouseIM(ϵ_d=μ) : AndersonIM(U=U, ϵ_d=μ)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		adt = reweighting!(lattice, mpsK, flat, mpsI, trunc=trunc)
		for band in 1:bands
			adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
		end

		cache = environments(lattice, adt)
		# 显式循环而非 comprehension（Julia 1.10 codegen bug，见 continuousbath/realtime.jl 注释）
		gt = Vector{ComplexF64}(undef, Nt+1)
		lt = Vector{ComplexF64}(undef, Nt+1)
		gτ = Vector{ComplexF64}(undef, Nτ+1)
		for k in 1:Nt+1
			gt[k] = -im * cached_gf(lattice, (GrassmannContourIndex(k, conj=false, branch=:+, band=1), GrassmannContourIndex(1, conj=true, branch=:+, band=1)), adt; cache=cache)
			lt[k] = im * cached_gf(lattice, (GrassmannContourIndex(1, conj=true, branch=:-, band=1), GrassmannContourIndex(k, conj=false, branch=:+, band=1)), adt; cache=cache)
		end
		for k in 1:Nτ+1
			gτ[k] = cached_gf(lattice, (GrassmannContourIndex(k, conj=false, branch=:τ, band=1), GrassmannContourIndex(1, conj=true, branch=:τ, band=1)), adt; cache=cache)
		end

		gt_ed = [independentbosons_greater_μ(spec, tj, β=β, μ=μ, U=U, bands=bands) for tj in ts]
		lt_ed = [independentbosons_lesser_μ(spec, tj, β=β, μ=μ, U=U, bands=bands) for tj in ts]
		gτ_ed = independentbosons_Gτ_μ(spec, β=β, μ=μ, Nτ=Nτ, U=U, bands=bands)

		@test relerr(gt, gt_ed) < 5.0e-2
		@test relerr(lt, lt_ed) < 5.0e-2
		@test relerr(gτ, gτ_ed) < 1.0e-2
	end
end
