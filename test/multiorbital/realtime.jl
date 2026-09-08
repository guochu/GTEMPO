@testset "Multiorbital, real time: Kanamori (norb = 2)" begin
	U = 1.0; J = 0.2; μ = -U / 2; ω = 1.0; α = 0.5
	β = 1.0; δt = 0.05; Nt = 5
	ts = 0:δt:(Nt*δt)
	# the 4-band real-time dynamics requires a large bond dimension
	trunc = truncdimcutoff(D=800, ϵ=1.0e-10)

	# ED reference with a separable thermal initial state
	H, a, adag, H0 = kanamori_ed(U=U, J=J, norb=2, μ=μ, ω=ω, α=α)
	gt_ed, lt_ed = greater_lesser_ed(H, a, adag, H0, ts, β)

	bath = fermionicbath(DiracDelta(ω=ω, α=α), β=β)
	model = KanamoriIM(U=U, J=J, norb=2, μ=μ)
	lat = GrassmannLattice(N=Nt, δt=δt, contour=:real, bands=4)
	lattice1 = similar(lat, bands=1)
	corr = correlationfunction(bath, lattice1)
	mpsI = hybriddynamics(lattice1, corr, trunc=trunc)
	Is = [fillband(lat, mpsI, band=b) for b in 1:4]

	mpsK = sysdynamics(lat, model, trunc=trunc)
	for band in 1:4
		mpsK = boundarycondition!(mpsK, lat, band=band, trunc=trunc)
	end
	mpsK = systhermalstate!(mpsK, lat, model, trunc=trunc, β=β)

	cache = environments(lat, mpsK, Is...)
	gt = [-im * cached_greater(lat, k, mpsK, Is..., band=1, cache=cache) for k in 1:lat.k]
	lt = [-im * cached_lesser(lat, k, mpsK, Is..., band=1, cache=cache) for k in 1:lat.k]
	@test relerr(gt, gt_ed) < 3.0e-2
	@test relerr(lt, lt_ed) < 3.0e-2
end
