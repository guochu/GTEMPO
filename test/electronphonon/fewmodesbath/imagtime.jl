@testset "Electron-phonon, imaginary time: independent bosons (phonon only)" begin
	μ = 0.5
	δτ = 0.1; N = 10; β = N * δτ
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)
	spec = DiracDelta(ω=1, α=0.5)

	for (U, bands) in ((0.0, 1), (1.0, 2))
		lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=bands)
		flat = FockLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=bands)
		bath = bosonicbath(spec, β=β)
		corr = correlationfunction(bath, flat)

		# fast and naive IF constructions agree
		mpsI = hybriddynamics(flat, corr, trunc=trunc)
		mpsI′ = hybriddynamics_naive(flat, corr, trunc=trunc)
		@test distance(mpsI, mpsI′) / norm(mpsI) < 1.0e-5

		model = (U == 0) ? ToulouseIM(μ=μ) : AndersonIM(U=U, μ=μ)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		adt = reweighting!(lattice, mpsK, flat, mpsI, trunc=trunc)
		for band in 1:bands
			adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
		end
		cache = environments(lattice, adt)
		g1 = cached_gf_fast(lattice, adt; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
		g2 = independentbosons_Gτ_μ(spec, β=β, μ=μ, Nτ=N, U=U, bands=bands)
		@test relerr(g1, g2) < 1.0e-2
	end

	# phonon-only IF built on the Grassmann lattice (retarded-interaction path)
	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, lattice)
	mpsI = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
	model = ToulouseIM(μ=μ)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	cache = environments(lattice, mpsK, mpsI)
	g1 = cached_gf_fast(lattice, mpsK, mpsI; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
	g2 = independentbosons_Gτ_μ(spec, β=β, μ=μ, Nτ=N)
	@test relerr(g1, g2) < 1.0e-2
end

@testset "Electron-phonon, imaginary time: fermionic + phonon bath" begin
	μ = 0.5
	δτ = 0.1; N = 10; β = N * δτ
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)
	ωf, αf = 1.0, 0.5
	ω₀, α = 1.0, 0.5

	for (U, bands) in ((0.0, 1), (1.0, 2))
		# ED reference for the combined model
		H, a, adag, H0 = mixedbath_ed(μ=μ, U=U, ωf=ωf, αf=αf, ω₀=ω₀, α=α, d=8)
		g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

		lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=bands)
		flat = FockLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=bands)

		# fermionic bath IF (Grassmann path, replicated per band)
		fbath = fermionicbath(DiracDelta(ω=ωf, α=αf), β=β)
		lattice1 = similar(lattice, bands=1)
		fcorr = correlationfunction(fbath, lattice1)
		mpsI_e = hybriddynamics(lattice1, fcorr, trunc=trunc)
		Is = [fillband(lattice, mpsI_e, band=b) for b in 1:bands]

		# phonon IF (Fock-lattice path, reweighted onto the Grassmann lattice)
		pbath = bosonicbath(DiracDelta(ω=ω₀, α=α), β=β)
		pcorr = correlationfunction(pbath, flat)
		mpsI_p = hybriddynamics(flat, pcorr, trunc=trunc)

		model = (U == 0) ? ToulouseIM(μ=μ) : AndersonIM(U=U, μ=μ)
		mpsK = sysdynamics(lattice, model, trunc=trunc)
		adt = reweighting!(lattice, mpsK, flat, mpsI_p, trunc=trunc)
		for band in 1:bands
			adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
		end

		cache = environments(lattice, adt, Is...)
		g = cached_gf_fast(lattice, adt, Is...; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
		# d = 8 phonon truncation in the ED reference limits the accuracy
		@test relerr(g, g_ed) < 2.0e-2
	end
end
