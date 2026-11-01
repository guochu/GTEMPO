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
	# the occupation observable is cell-centred (its density pair spans the
	# (j+1, j) legs), so node k corresponds to t = (k - ½)δt: sampling the ED
	# reference at i·δt misaligns it by half a step and leaves a spurious
	# O(δt) error (0.23 at δt=0.05; with the half shift 0.014 → 0.004, i.e.
	# second-order convergent). The current is a cell rate and keeps the
	# end-of-cell time.
	for i in 1:N
		ρ = timeevo(ρ₀, h, -im*(i - 0.5)*δt, cache_ed)
		push!(ns_ed, real(ρ[1, 1]))
		ρ_end = timeevo(ρ₀, h, -im*i*δt, cache_ed)
		push!(currents_ed, sum(jobs .* ρ_end))
	end

	exact_model = AndersonIM(ϵ_d=ϵ_d, U=0)
	lattice_o = GrassmannLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lattice_o)

	lattice = similar(lattice_o, N=0)
	mpsI = vacuumstate(lattice)
	mpsK = vacuumstate(lattice)
	for k in 2:N+1
		lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
		@test k == timesteps(mpsI, lattice) == lattice.k
		mpsI = hybriddynamicsstepper!(mpsI, lattice, corr, trunc=trunc)
		mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=trunc)
	end
	# measure on the final lattice with the boundary condition applied once:
	# applying it to the growing (future-truncated) lattice at every step
	# conflates the online evolution with the truncated-contour BC, which the
	# vertex convention no longer tolerates (0.11 vs 0.014 here)
	mpsI′ = boundarycondition(mpsI, lattice)
	ns = [real(occupation(lattice, n, mpsK, mpsI′, branch=:+)) for n in 1:N]
	currents = [electriccurrent_fast(lattice, corr, k, mpsK, mpsI′) for k in 2:N+1]
	@test length(ns) == length(currents) == N
	@test norm(ns - ns_ed) / norm(ns_ed) < rtol
	@test norm(currents - currents_ed) / norm(currents_ed) < rtol
end
