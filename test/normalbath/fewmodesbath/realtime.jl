@testset "Few-mode bath, real time: orderings & algorithms" begin
	μ = 0.7; ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 8
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=60, ϵ=1.0e-10)

	# ED reference: one fermionic bath mode, separable thermal initial state
	H, a, adag, H0 = singlemode_ed(μ=μ, U=0, bathspecs=[(ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

	spec = DiracDelta(ω=ω, α=α)
	bath = fermionicbath(spec, β=β)
	model = AndersonIM(U=0, μ=μ)
	for ordering in real_orderings
		lat = GrassmannLattice(N=Nt, δt=δt, contour=:real, ordering=ordering)
		corr = correlationfunction(bath, lat)
		# all algorithms on the default ordering, the default algorithm on all orderings
		algs = (ordering == real_orderings[1]) ? if_algs(trunc) :
			   [("PartialIF", PartialIF(trunc=trunc))]
		for (name, alg) in algs
			mpsI = hybriddynamics(lat, corr, alg)
			mpsK = sysdynamics(lat, model, trunc=trunc)
			mpsK = boundarycondition!(mpsK, lat)
			mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)
			gt, lt = gtlt_series(lat, mpsK, mpsI)
			@test relerr(gt, gt_ed) < 3.0e-2
			@test relerr(lt, lt_ed) < 3.0e-2
		end
	end
end

@testset "Few-mode bath, real time: two bands with U" begin
	μ = 0.7; U = 1.0
	β = 1.0; δt = 0.05; Nt = 10
	ts = 0:δt:(Nt*δt)
	trunc = truncdimcutoff(D=80, ϵ=1.0e-10)
	ω, α = 1.0, 0.5

	# two bands with U: ED reference (16-dimensional ED), one bath mode per band
	H, a, adag, H0 = singlemode_ed(μ=μ, U=U, bathspecs=[(ω, α), (ω, α)])
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real, bands=2)
	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	mpsK, Is = fermionic_setup(lat, bath, AndersonIM(U=U, μ=μ), trunc, β=β)
	gt, lt = gtlt_series(lat, mpsK, Is...)
	@test relerr(gt, gt_ed) < 3.0e-2
	@test relerr(lt, lt_ed) < 3.0e-2
end
