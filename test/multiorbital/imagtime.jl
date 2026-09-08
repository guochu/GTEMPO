@testset "Multiorbital, imaginary time: two bands with U" begin
	μ = 0.7; U = 1.0; ω = 1.0; α = 0.5
	δτ = 0.1; N = 8; β = N * δτ
	# the band-interleaved ordering (last of imag_orderings) needs a larger bond dimension
	trunc = truncdimcutoff(D=160, ϵ=1.0e-10)

	H, a, adag, H0 = singlemode_ed(μ=μ, U=U, bathspecs=[(ω, α), (ω, α)])
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = AndersonIM(U=U, μ=μ)
	for ordering in imag_orderings
		lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag, ordering=ordering, bands=2)
		mpsK, Is = fermionic_setup(lat, bath, model, trunc)
		g = gτ_series(lat, mpsK, Is...)
		@test relerr(g, g_ed) < 1.0e-2
	end
end

@testset "Multiorbital, imaginary time: Kanamori (norb = 2)" begin
	U = 1.0; J = 0.2; μ = -U / 2; ω = 1.0; α = 0.5
	δτ = 0.1; N = 6; β = N * δτ
	# the 4-band interacting dynamics requires a large bond dimension
	trunc = truncdimcutoff(D=400, ϵ=1.0e-10)

	# ED reference: 4 impurity bands + 4 bath modes (256-dimensional space)
	H, a, adag, H0 = kanamori_ed(U=U, J=J, norb=2, μ=μ, ω=ω, α=α)
	g_ed = gτ_ed(H, a, adag, 0:δτ:β, β)

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = KanamoriIM(U=U, J=J, norb=2, μ=μ)
	lat = GrassmannLattice(N=N, δτ=δτ, contour=:imag, bands=4)
	lattice1 = similar(lat, bands=1)
	corr = correlationfunction(bath, lattice1)
	mpsI = hybriddynamics(lattice1, corr, trunc=trunc)
	Is = [fillband(lat, mpsI, band=b) for b in 1:4]

	mpsK = sysdynamics(lat, model, trunc=trunc)
	for band in 1:4
		mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
	end
	cache = environments(lat, mpsK, Is...)
	g = cached_gf_fast(lat, mpsK, Is...; c1=false, c2=true, b1=:τ, b2=:τ, cache=cache)
	@test relerr(g, g_ed) < 1.0e-2

	# band 2 Green's function as an additional check
	g2 = cached_gf_fast(lat, mpsK, Is...; c1=false, c2=true, b1=:τ, b2=:τ, band=2, cache=cache)
	cs = fermion_annihilators(8)
	g2_ed = gτ_ed(H, cs[2], cs[2]', 0:δτ:β, β)
	@test relerr(g2, g2_ed) < 1.0e-2
end
