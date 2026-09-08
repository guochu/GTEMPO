# real-time particle currents on the Keldysh contour, benchmarked against the
# exact-diagonalization reference solutions of ImpurityModelBase:
# - the (single-bath) Toulouse model with a discretized continuous bath
# - the boundary-driven double-bath resonant level

@testset "Continuous bath, real time: particle current (Toulouse ED)" begin
	β = 1.0
	ϵ_d = -1.0
	N = 8
	δt = 0.03
	dw = 0.05
	trunc = truncdimcutoff(D=150, ϵ=1.0e-9)
	rtol = 5.0e-2

	for μ in (0.0, 1.5)
		# ED reference: discretize the continuous bath (Toulouse model)
		bath = fermionicbath(semicircular(5.0), β=β, μ=μ)
		model = Toulouse(discretebath(bath, δw=dw), ϵ_d=ϵ_d)
		ρ₀ = separablecdm(model, 0)
		h = cmatrix(model)
		cache_ed = eigencache(transpose(h))
		observer = particlecurrent_cmatrix(model)
		currents_ed = [sum(observer .* timeevo(ρ₀, h, -im*i*δt, cache_ed)) for i in 1:N]

		# GTEMPO: plain series
		lattice = GrassmannLattice(N=N, δt=δt, contour=:real)
		corr = correlationfunction(bath, lattice)
		mpsI = hybriddynamics(lattice, corr, trunc=trunc)
		mpsI = boundarycondition!(mpsI, lattice)
		mpsK = sysdynamics(lattice, AndersonIM(μ=ϵ_d, U=0), trunc=trunc)

		currents = electriccurrent(lattice, corr, mpsK, mpsI)
		@test length(currents) == N
		@test norm(currents - currents_ed) / norm(currents_ed) < rtol

		# MPO-based variant at the last time step
		curr_fast = electriccurrent_fast(lattice, corr, N + 1, mpsK, mpsI)
		@test abs(curr_fast - currents_ed[end]) / norm(currents_ed) < rtol

		# cached evaluation agrees with the plain one
		cache = environments(lattice, mpsK, mpsI)
		currents_cached = cached_electriccurrent(lattice, corr, mpsK, mpsI, cache=cache)
		@test norm(currents_cached - currents) / norm(currents) < 1.0e-7
	end
end

@testset "Continuous bath, real time: boundary-driven double bath (BoundaryDriving ED)" begin
	β = 2.0
	V = 1.0
	ϵ_d = -0.5
	N = 8
	δt = 0.05
	dw = 0.05
	trunc = truncdimcutoff(D=100, ϵ=1.0e-9)
	rtol = 6.0e-2

	leftbath = fermionicbath(semicircular(5.0), β=β, μ=V / 2)
	rightbath = fermionicbath(semicircular(5.0), β=β, μ=-V / 2)

	# ED reference: single-site system driven by the two discretized baths
	model = BoundaryDriving([ϵ_d;;], discretebath(leftbath, δw=dw), discretebath(rightbath, δw=dw))
	ρ₀ = separablecdm(model, zeros(1, 1))
	h = cmatrix(model)
	cache_ed = eigencache(transpose(h))
	jl = leftparticlecurrent_cmatrix(model)
	jr = rightparticlecurrent_cmatrix(model)
	currents_l = [sum(jl .* timeevo(ρ₀, h, -im*i*δt, cache_ed)) for i in 1:N]
	currents_r = [sum(jr .* timeevo(ρ₀, h, -im*i*δt, cache_ed)) for i in 1:N]

	# GTEMPO: the total IF is the sum of the two lead IFs on the same impurity band;
	# the current through each lead uses its own correlation function
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real)
	lcorr = correlationfunction(leftbath, lattice)
	rcorr = correlationfunction(rightbath, lattice)
	corr = lcorr + rcorr
	mpsI = hybriddynamics(lattice, corr, trunc=trunc)
	mpsI = boundarycondition!(mpsI, lattice)
	mpsK = sysdynamics(lattice, AndersonIM(μ=ϵ_d, U=0), trunc=trunc)

	cl = electriccurrent(lattice, lcorr, mpsK, mpsI)
	cr = electriccurrent(lattice, rcorr, mpsK, mpsI)
	@test norm(cl - currents_l) / norm(currents_l) < rtol
	@test norm(cr - currents_r) / norm(currents_r) < rtol

	# cached MPO variant agrees with the plain series
	cache = environments(lattice, mpsK, mpsI)
	cl2 = cached_electriccurrent_fast(lattice, lcorr, mpsK, mpsI, cache=cache)
	cr2 = cached_electriccurrent_fast(lattice, rcorr, mpsK, mpsI, cache=cache)
	@test norm(cl2 - cl) / norm(cl) < 1.0e-6
	@test norm(cr2 - cr) / norm(cr) < 1.0e-6
end
