println("------------------------------------")
println("|          BCS Real-time           |")
println("------------------------------------")




@testset "BCS Real time evolution" begin
	tol = 1.0e-5
	rtol = 1.0e-2
	δt=0.05
	Nt = 10
	β = 10
	t = Nt * δt
	chi = 60
	chi2 = 200

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10)

	U = 1
	# ϵ_d = -0.7
	ϵ_d = 0.8
	ω = 1
	α = 0.5

	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:Keldysh, bands=2)
	model = AndersonIM(U=U, μ=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end

	# trivial case
	Δ = 0
	bath = fermionicbath(semicircular(t=1), β=β)
	corr = correlationfunction(bath, lattice)
	mpsI1 = hybriddynamics(lattice, corr, band=1, trunc=trunc)
	mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)

	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	g1 = [-im*cached_greater(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.k]
	g2 = [-im*cached_lesser(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.k]

	bath2 = bcsbath(bath)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc2)
	cache = environments(lattice, mpsK, mpsI)
	g1′ = [-im*cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	g2′ = [-im*cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	@test norm(g1-g1′) / norm(g1) < rtol
	@test norm(g2-g2′) / norm(g2) < 5*rtol

	mpsK = systhermalstate!(mpsK, lattice, model, trunc=trunc, β=β)

	# compare with ED
	Δ = 0.6
	bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc2)
	mpsI′ = hybriddynamics!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) < tol
	cache = environments(lattice, mpsK, mpsI)
	g1 = [-im*cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	g2 = [-im*cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]


	H, a, adag, H0 = bcs_operators(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)

	ρ = exp(-β*H0)
	cache = eigencache(H)
	g1′ = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, cache, reverse = false)
	g2′ = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, cache, reverse = true)

	@test norm(g1-g1′) / norm(g1) < rtol
	@test norm(g2-g2′) / norm(g2) < rtol

	# complex gap
	Δ = 0.3 + 0.4*im
	bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc2)
	mpsI′ = hybriddynamics!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) < tol
	cache = environments(lattice, mpsK, mpsI)
	g1 = [-im*cached_greater(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]
	g2 = [-im*cached_lesser(lattice, i, mpsK, mpsI, cache=cache) for i in 1:lattice.k]


	H, a, adag, H0 = bcs_operators(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)

	ρ = exp(-β*H0)
	cache = eigencache(H)
	g1′ = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, cache, reverse = false)
	g2′ = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, cache, reverse = true)

	@test norm(g1-g1′) / norm(g1) < rtol
	@test norm(g2-g2′) / norm(g2) < rtol
end
