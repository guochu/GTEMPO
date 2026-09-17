@testset "API: hybriddynamics (imaginary time)" begin
	β = 0.3
	N = 3
	trunc = truncdimcutoff(D=30, ϵ=1.0e-10)
	bath = fermionicbath(DiracDelta(ω=1, α=0.5), β=β)

	for ordering in imag_orderings
		lat = GrassmannLattice(N=N, δτ=β/N, contour=:imag, ordering=ordering)
		corr = correlationfunction(bath, lat)
		# default (PartialIF) construction via trunc keyword
		mpsI = hybriddynamics(lat, corr, trunc=trunc)
		@test mpsI isa GrassmannMPS
		@test length(mpsI) == length(lat)
		# explicit PartialIF
		mpsI1 = hybriddynamics(lat, corr, PartialIF(trunc=trunc))
		# naive construction
		mpsI2 = hybriddynamics_naive(lat, corr, trunc=trunc)
		# XTRGIF, ExactTTIIF, TDVPIF
		mpsI3 = hybriddynamics(lat, corr, XTRGIF(k=2, algmult=SVDCompression(trunc), verbosity=0))
		mpsI4 = hybriddynamics(lat, corr, ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0))
		mpsI5 = hybriddynamics(lat, corr, TDVPIF(trunc=trunc, δ=0.1, verbosity=0))
		for mps in (mpsI1, mpsI2, mpsI3, mpsI4, mpsI5)
			@test norm(mps) > 0
			@test relerr(mps, mpsI) < 5.0e-2
		end
		# in-place variants
		mpsI6 = hybriddynamics!(vacuumstate(lat), lat, corr, trunc=trunc)
		@test relerr(mpsI6, mpsI) < 1.0e-6
	end

	# in-place ExactTTIIF (default ordering): the IF is merged into the state
	# passed in, e.g. the impurity dynamics obtained from sysdynamics
	lat = GrassmannLattice(N=N, δτ=β/N, contour=:imag)
	corr = correlationfunction(bath, lat)
	alg = ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0)
	mpsI4 = hybriddynamics(lat, corr, alg)
	@test relerr(hybriddynamics!(vacuumstate(lat), lat, corr, alg), mpsI4) < 1.0e-12
	model = ToulouseIM(ϵ_d=0.5)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = hybriddynamics!(mpsK, lat, corr, alg)
	mpsK_ref = mult(sysdynamics(lat, model, trunc=trunc), mpsI4, trunc=trunc)
	@test relerr(mpsK, mpsK_ref) < 1.0e-6

	# multi-band: construct on a single-band lattice, then fillband
	lat = GrassmannLattice(N=N, δτ=β/N, contour=:imag, bands=2)
	lat1 = similar(lat, bands=1)
	corr = correlationfunction(bath, lat1)
	mpsI = hybriddynamics(lat1, corr, trunc=trunc, band=1)
	mpsI1 = fillband(lat, mpsI, band=1)
	mpsI2 = fillband(lat, mpsI, band=2)
	@test length(mpsI1) == length(lat)
	@test relerr(mpsI2, swapband(mpsI1, lat, 1, 2, trunc=trunc)) < 1.0e-8

	# ExactTTIIF builds its terms on a single-band lattice internally and
	# expands them to the full lattice via fillband
	mpsI3 = hybriddynamics(lat, corr, ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0))
	mpsI3_ref = fillband(lat, hybriddynamics(lat1, corr, ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0)), band=1)
	@test relerr(mpsI3, mpsI3_ref) < 1.0e-12
end

@testset "API: hybriddynamics (real time)" begin
	δt = 0.05
	N = 3
	trunc = truncdimcutoff(D=30, ϵ=1.0e-10)
	bath = fermionicbath(DiracDelta(ω=1, α=0.5), β=1.0)

	for ordering in real_orderings
		lat = GrassmannLattice(N=N, δt=δt, contour=:real, ordering=ordering)
		corr = correlationfunction(bath, lat)
		mpsI = hybriddynamics(lat, corr, trunc=trunc)
		@test mpsI isa GrassmannMPS
		@test length(mpsI) == length(lat)
		mpsI1 = hybriddynamics(lat, corr, PartialIF(trunc=trunc))
		mpsI2 = hybriddynamics_naive(lat, corr, trunc=trunc)
		mpsI3 = hybriddynamics(lat, corr, XTRGIF(k=2, algmult=SVDCompression(trunc), verbosity=0))
		mpsI4 = hybriddynamics(lat, corr, ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0))
		mpsI5 = hybriddynamics(lat, corr, TDVPIF(trunc=trunc, δ=0.1, verbosity=0))
		for mps in (mpsI1, mpsI2, mpsI3, mpsI4, mpsI5)
			@test norm(mps) > 0
			@test relerr(mps, mpsI) < 5.0e-2
		end
	end

	# in-place ExactTTIIF (default ordering): merge the IF into a K obtained
	# from sysdynamics, on the real-time (Keldysh) contour
	lat = GrassmannLattice(N=N, δt=δt, contour=:real)
	corr = correlationfunction(bath, lat)
	alg = ExactTTIIF(algmult=SVDCompression(trunc), verbosity=0)
	mpsI4 = hybriddynamics(lat, corr, alg)
	@test relerr(hybriddynamics!(vacuumstate(lat), lat, corr, alg), mpsI4) < 1.0e-12
	model = ToulouseIM(ϵ_d=0.5)
	mpsK = sysdynamics(lat, model, trunc=trunc)
	mpsK = hybriddynamics!(mpsK, lat, corr, alg)
	mpsK_ref = mult(sysdynamics(lat, model, trunc=trunc), mpsI4, trunc=trunc)
	@test relerr(mpsK, mpsK_ref) < 1.0e-3
end

@testset "API: hybriddynamics (mixed time)" begin
	β = 0.3
	trunc = truncdimcutoff(D=30, ϵ=1.0e-10)
	bath = fermionicbath(DiracDelta(ω=1, α=0.5), β=β)

	for ordering in mixed_orderings
		lat = GrassmannLattice(Nt=3, δt=0.05, Nτ=3, δτ=0.1, contour=:mixed, ordering=ordering)
		corr = correlationfunction(bath, lat)
		mpsI = hybriddynamics(lat, corr, trunc=trunc)
		@test mpsI isa GrassmannMPS
		@test length(mpsI) == length(lat)
		mpsI2 = hybriddynamics_naive(lat, corr, trunc=trunc)
		@test norm(mpsI) > 0
		@test relerr(mpsI, mpsI2) < 5.0e-2
		# per-band construction
		mpsI3 = hybriddynamics(lat, corr, trunc=trunc, band=1)
		@test relerr(mpsI3, mpsI) < 1.0e-8
	end
end

@testset "API: influence functional internals" begin
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	bath = fermionicbath(DiracDelta(ω=1, α=0.5), β=0.4)

	# algorithm type hierarchy
	@test PartialIF() isa InfluenceFunctionalAlgorithm
	@test XTRGIF() isa InfluenceFunctionalAlgorithm
	@test ExactTTIIF() isa InfluenceFunctionalAlgorithm
	@test TDVPIF() isa InfluenceFunctionalAlgorithm

	# imaginary time: influenceoperators / influenceoperatorsteppers / influenceoperatorstepper
	lat = GrassmannLattice(N=3, δτ=0.1, contour=:imag)
	corr = correlationfunction(bath, lat)
	ops = influenceoperators(lat, corr)
	@test length(ops) == 1 && ops[1] isa MPO && length(ops[1]) == length(lat)
	ops1 = influenceoperatorsteppers(lat, corr, 0.05, WII(tol=1.0e-14, maxiter=100000))
	@test length(ops1) == 1 && ops1[1] isa MPO
	ops2 = influenceoperatorsteppers(lat, corr, 0.05, ComplexStepper(WII(tol=1.0e-14, maxiter=100000)))
	@test length(ops2) == 2 && all(op -> op isa MPO, ops2)
	mps = influenceoperatorstepper(lat, corr, 0.05, WII(tol=1.0e-14, maxiter=100000), SVDCompression(trunc))
	@test mps isa GrassmannMPS && length(mps) == length(lat)
	@test relerr(mps, ops1[1] * vacuumstate(lat)) < 1.0e-8

	# partialif_hybrid (analytic MPO) vs partialif_hybrid_naive (GTerm exponentials)
	k = lat.k - 1
	mps1 = vacuumstate(lat)
	for i in 1:k
		mps1 = mult!(mps1, partialif_hybrid(lat, i + 1, [0; view(corr.data, i, 1:k)]), trunc=trunc)
	end
	@test relerr(mps1, hybriddynamics(lat, corr, trunc=trunc)) < 1.0e-8
	a = partialif_hybrid(lat, 2, [0; view(corr.data, 1, 1:k)])
	b = partialif_hybrid_naive(lat, 2, [0; view(corr.data, 1, 1:k)], trunc=trunc)
	@test relerr(a, b) < 1.0e-8
	# hybriddynamics_naive! (in place) agrees with hybriddynamics_naive
	@test relerr(hybriddynamics_naive!(vacuumstate(lat), lat, corr, trunc=trunc),
				 hybriddynamics_naive(lat, corr, trunc=trunc)) < 1.0e-12

	# real time: stepwise evolution via hybriddynamicsstepper!
	lat_r = GrassmannLattice(N=3, δt=0.05, contour=:real)
	corr_r = correlationfunction(bath, lat_r)
	lattice = similar(lat_r, N=0)
	mpsI = vacuumstate(lattice)
	while lattice.N < lat_r.N
		lattice, mpsI = makestep(lattice, mpsI)
		mpsI = hybriddynamicsstepper!(mpsI, lattice, corr_r, trunc=trunc)
	end
	@test timesteps(mpsI, lattice) == lattice.k
	# copy and in-place steppers agree
	@test distance(hybriddynamicsstepper(mpsI, lattice, corr_r, trunc=trunc),
				   hybriddynamicsstepper!(copy(mpsI), lattice, corr_r, trunc=trunc)) == 0
	# the stepped IF equals the static construction on the full lattice
	# (small accumulated truncation differences between the two paths)
	@test relerr(mpsI, hybriddynamics(lat_r, corr_r, trunc=trunc)) < 1.0e-5

	# second-order lattice: finalized IF built on a copy of the carried state
	lat_r2 = GrassmannLattice(N=2, δt=0.05, contour=:real, order=2)
	corr_r2 = correlationfunction(bath, lat_r2)
	lattice = similar(lat_r2, N=0)
	mpsI2 = vacuumstate(lattice)
	for k2 in 2:lat_r2.N+1
		lattice, mpsI2 = makestep(lattice, mpsI2)
		mpsI2 = hybriddynamicsstepper!(mpsI2, lattice, corr_r2, finalize=false, trunc=trunc)
	end
	mpsI2f = hybriddynamicsstepper!(copy(mpsI2), lattice, corr_r2, finalize=true, trunc=trunc)
	@test length(mpsI2f) == length(lattice) && norm(mpsI2f) > 0

	# retarded interaction (bosonic bath): naive construction, copy vs in-place
	bblat = GrassmannLattice(N=3, δτ=0.1, contour=:imag)
	bbath = bosonicbath(DiracDelta(ω=1, α=0.5), β=1.0)
	bcorr = correlationfunction(bbath, bblat)
	V1 = retardedinteractdynamics_naive(bblat, bcorr, trunc=trunc)
	V2 = retardedinteractdynamics_naive!(vacuumstate(bblat), bblat, bcorr, trunc=trunc)
	@test length(V1) == length(bblat) && norm(V1) > 0
	@test distance(V1, V2) < 1.0e-10
end
