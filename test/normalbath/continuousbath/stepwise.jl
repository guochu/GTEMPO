# stepwise (online) evolution on the real-time lattice: makestep + hybriddynamicsstepper!
# + boundarycondition + sysdynamicsstepper!, benchmarked against the Toulouse ED reference
# of ImpurityModelBase
#
# NOTE: the second-order stepwise recipe follows the archived test: the finalized
# (full 2k-2 row) IF is built on a copy of the carried-over MPS for the observables,
# then the carried-over MPS itself receives the full two-row (2k-2 + 2k-1) increment.
# Calling hybriddynamicsstepper! in-place on the same MPS for both calls would apply
# the 2k-2 row twice.

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

	exact_model = AndersonIM(μ=ϵ_d, U=0)
	for order in (1, 2)
		lattice_o = GrassmannLattice(N=N, δt=δt, contour=:real, order=order)
		corr = correlationfunction(bath, lattice_o)

		lattice = similar(lattice_o, N=0)
		mpsI = vacuumstate(lattice)
		mpsK = vacuumstate(lattice)
		ns = Float64[]
		currents = ComplexF64[]
		for k in 2:N+1
			lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
			@test k == timesteps(mpsI, lattice) == lattice.k
			if order == 1
				mpsI = hybriddynamicsstepper!(mpsI, lattice, corr, trunc=trunc)
				mpsI′ = boundarycondition(mpsI, lattice)
				mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=trunc)
				push!(ns, real(occupation(lattice, k - 1, mpsK, mpsI′, branch=:+)))
				push!(currents, electriccurrent_fast(lattice, corr, k, mpsK, mpsI′))
			else
				# see the NOTE above: the finalized IF goes to a copy of the carried MPS
				mpsI2 = hybriddynamicsstepper!(copy(mpsI), lattice, corr, finalize=true, trunc=trunc)
				mpsI2 = boundarycondition(mpsI2, lattice)
				mpsI = hybriddynamicsstepper!(mpsI, lattice, corr, finalize=false, trunc=trunc)
				mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=trunc)
				cache = environments(lattice, mpsK, mpsI2)
				# for the second-order lattice the cached observables are defined at the
				# current (largest) time step only
				push!(ns, real(cached_occupation(lattice, mpsK, mpsI2, cache=cache)))
				push!(currents, cached_electriccurrent_fast(lattice, corr, mpsK, mpsI2, cache=cache))
			end
		end
		@test length(ns) == length(currents) == N
		@test norm(ns - ns_ed) / norm(ns_ed) < rtol
		@test norm(currents - currents_ed) / norm(currents_ed) < rtol
	end
end
