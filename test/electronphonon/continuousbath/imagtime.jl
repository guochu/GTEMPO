@testset "Electron-phonon (continuous bath), imaginary time: independent bosons" begin
	μ = 0.5
	δτ = 0.1; N = 10; β = N * δτ
	trunc = truncdimcutoff(D=100, ϵ=1.0e-10)
	spec = Leggett(d=3, ωc=1)

	for (U, bands) in ((0.0, 1), (1.0, 2))
		lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=bands)
		flat = FockLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=bands)
		bath = bosonicbath(spec, β=β)
		corr = correlationfunction(bath, flat)

		mpsI = hybriddynamics(flat, corr, trunc=trunc)
		mpsI′ = hybriddynamics_naive(flat, corr, trunc=trunc)
		@test distance(mpsI, mpsI′) / norm(mpsI) < 1.0e-5

		model = AndersonIM(U=U, μ=μ)
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
end
