@testset "Continuous bath, real time: Toulouse (algorithms)" begin
	μ = 1.25 * pi
	spec = spectrum_func()
	β = 1.0; δt = 0.03; Nt = 6
	ts = [i*δt for i in 0:Nt]
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)

	bath = fermionicbath(spec, β=β, μ=0)
	# equilibrium Toulouse Green's function (GTEMPO convention: im * toulouse_Gt)
	gt_ref = [im * toulouse_Gt(bath, tj; ϵ_d=μ) for tj in ts]

	model = AndersonIM(U=0, μ=μ)
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real)
	corr = correlationfunction(bath, lat)
	# TDVPIF (real time, continuous bath) is the most expensive construction;
	# the full algorithm set is covered by the single-mode bath tests
	for (name, alg) in [("PartialIF", PartialIF(trunc=trunc)),
					    ("ExactTTIIF", ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0))]
		mpsI = hybriddynamics(lat, corr, alg)
		mpsK = sysdynamics(lat, model, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lat)
		mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)
		gt = [gf(lat, (ContourIndex(k, conj=false, branch=:+, band=1), ContourIndex(1, conj=true, branch=:+, band=1)), mpsK, mpsI; Z=integrate(lat, mpsK, mpsI)) for k in 1:lat.k]
		@test relerr(gt, gt_ref) < 5.0e-2
	end
end
