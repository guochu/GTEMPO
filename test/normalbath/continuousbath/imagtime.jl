@testset "Continuous bath, imaginary time: Toulouse (algorithms)" begin
	μ = 1.25 * pi
	spec = spectrum_func()
	δτ = 0.03; N = 10; β = N * δτ
	τs = collect(0:δτ:β)
	trunc = truncdimcutoff(D=120, ϵ=1.0e-6)

	# two independent references: analytic Toulouse formula and the
	# discretized-bath Toulouse ED model
	bath = fermionicbath(spec, β=β, μ=0)
	g_ana = [toulouse_Gτ(bath, τ; ϵ_d=μ) for τ in τs]
	g_disc = toulouse_Gτ(Toulouse(discretebath(bath, δw=0.2), ϵ_d=μ), τs)
	@test relerr(g_disc, g_ana) < 1.0e-2

	model = AndersonIM(U=0, μ=μ)
	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	corr = correlationfunction(bath, lat)
	for (name, alg) in if_algs(trunc)
		mpsI = hybriddynamics(lat, corr, alg)
		mpsK = sysdynamics(lat, model, trunc=trunc)
		mpsK = boundarycondition!(mpsK, lat)
		g = gτ_series(lat, mpsK, mpsI)
		@test relerr(g, g_ana) < 1.0e-2
	end
end
