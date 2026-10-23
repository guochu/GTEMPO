# stepwise (online) evolution on the real-time lattice: makestep + hybriddynamicsstepper!
# + boundarycondition + sysdynamicsstepper!, benchmarked against the Toulouse ED reference
# of ImpurityModelBase

@testset "Continuous bath, real time: stepwise evolution (Toulouse ED)" begin
	β = 1.0
	μ = 1.0
	ϵ_d = -1.0
	N = 8
	δt = 0.05
	dw = 0.05
	trunc = truncdimcutoff(D=100, ϵ=1.0e-9)
	rtol = 5.0e-2

	bath = fermionicbath(semicircular(5.0), β=β, μ=μ)

	# ED reference
	model_ed = Toulouse(discretebath(bath, δw=dw), ϵ_d=ϵ_d)
	ρ₀ = separablecdm(model_ed, 0)
	h = cmatrix(model_ed)
	cache_ed = eigencache(transpose(h))
	jobs = particlecurrent_cmatrix(model_ed)
	ns_ed = Float64[]
	currents_ed = ComplexF64[]
	for i in 1:N
		ρ = timeevo(ρ₀, h, -im*i*δt, cache_ed)
		push!(ns_ed, real(ρ[1, 1]))
		push!(currents_ed, sum(jobs .* ρ))
	end

	exact_model = AndersonIM(ϵ_d=ϵ_d, U=0)
	lattice_o = GrassmannLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lattice_o)

	lattice = similar(lattice_o, N=0)
	mpsI = vacuumstate(lattice)
	mpsK = vacuumstate(lattice)
	ns = Float64[]
	currents = ComplexF64[]
	for k in 2:N+1
		lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
		@test k == timesteps(mpsI, lattice) == lattice.k
		mpsI = hybriddynamicsstepper!(mpsI, lattice, corr, trunc=trunc)
		mpsI′ = boundarycondition(mpsI, lattice)
		mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=trunc)
		push!(ns, real(occupation(lattice, k - 1, mpsK, mpsI′, branch=:+)))
		push!(currents, electriccurrent_fast(lattice, corr, k, mpsK, mpsI′))
	end
	@test length(ns) == length(currents) == N
	@test norm(ns - ns_ed) / norm(ns_ed) < rtol
	@test norm(currents - currents_ed) / norm(currents_ed) < rtol
end
