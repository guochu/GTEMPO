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

	# multi-band: construct on a single-band lattice, then fillband
	lat = GrassmannLattice(N=N, δτ=β/N, contour=:imag, bands=2)
	lat1 = similar(lat, bands=1)
	corr = correlationfunction(bath, lat1)
	mpsI = hybriddynamics(lat1, corr, trunc=trunc, band=1)
	mpsI1 = fillband(lat, mpsI, band=1)
	mpsI2 = fillband(lat, mpsI, band=2)
	@test length(mpsI1) == length(lat)
	@test relerr(mpsI2, swapband(mpsI1, lat, 1, 2, trunc=trunc)) < 1.0e-8
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
