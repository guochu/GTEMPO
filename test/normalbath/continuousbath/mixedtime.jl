@testset "Continuous bath, mixed time: Toulouse" begin
	μ = 1.25 * pi
	spec = spectrum_func()
	δτ = 0.02; Nτ = 15; β = Nτ * δτ
	δt = 0.02; Nt = 8
	ts = [i*δt for i in 0:Nt]
	τs = collect(0:δτ:β)
	trunc = truncdimcutoff(D=100, ϵ=1.0e-6)

	bath = fermionicbath(spec, β=β, μ=0)
	# Toulouse ED references on the discretized bath
	ed_model = Toulouse(discretebath(bath, δw=0.1), ϵ_d=μ)
	gt_ref, lt_ref = toulouse_greater_lesser(ed_model, ts)
	gt_ref, lt_ref = im * gt_ref, -im * lt_ref
	gτ_ref = toulouse_Gτ(ed_model, τs)

	model = ToulouseIM(μ=μ)
	lat = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed)
	corr = correlationfunction(bath, lat)
	mpsI = hybriddynamics(lat, corr, trunc=trunc)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = boundarycondition!(mpsK, lat)
	cache = environments(lat, mpsK, mpsI)
	gt = [cached_greater(lat, k, mpsK, mpsI, cache=cache) for k in 1:lat.kt]
	lt = [cached_lesser(lat, k, mpsK, mpsI, cache=cache) for k in 1:lat.kt]
	gτ = cached_gf_fast(lat, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
	gτ[end] = 1 - gτ[1]
	@test relerr(gt, gt_ref) < 2.0e-2
	@test relerr(lt, lt_ref) < 2.0e-2
	@test relerr(gτ, gτ_ref) < 2.0e-2
end
