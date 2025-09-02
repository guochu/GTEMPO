println("------------------------------------")
println("|        BCS Imaginary-time        |")
println("------------------------------------")




@testset "BCS Imaginary time evolution" begin
	tol = 1.0e-5
	rtol = 1.0e-2
	δτ=0.1
	Nτ = 10
	β = Nτ * δτ
	chi = 60

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	U = 1
	# ϵ_d = -0.7
	ϵ_d = 0.8
	ω = 1
	α = 0.5

	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
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
	g1 = cached_Gτ(lattice, mpsK, mpsI1, mpsI2, cache=cache)

	bath2 = bcsbath(bath)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics_naive!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc)
	cache = environments(lattice, mpsK, mpsI)
	g2 = cached_Gτ(lattice, mpsK, mpsI, cache=cache)
	@test norm(g1-g2) / norm(g1) < rtol


	# compare with ED
	Δ = 0.6
	bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc)
	mpsI′ = hybriddynamics!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) < tol
	cache = environments(lattice, mpsK, mpsI)
	g1 = cached_Gτ(lattice, mpsK, mpsI, cache=cache)

	H, a, adag, H0 = bcs_operators(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)

	g2 = correlation_2op_1τ(H, a, adag, 0:δτ:β, β=β)

	@test norm(g1-g2) / norm(g2) < rtol

	# complex gap
	Δ = 0.2 + 0.45*im
	mpsK = complex(mpsK)
	bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc)
	mpsI′ = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) < tol
	cache = environments(lattice, mpsK, mpsI)
	g1 = cached_Gτ(lattice, mpsK, mpsI, cache=cache)

	H, a, adag, H0 = bcs_operators(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)

	g2 = correlation_2op_1τ(H, a, adag, 0:δτ:β, β=β)

	@test norm(g1-g2) / norm(g2) < rtol

	# quadratic ED
	Δ = 0.2 + 0.45*im
	δτ=0.05
	β = Nτ * δτ
	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
	model = AndersonIM(U=0, μ=-ϵ_d)
	mpsK = sysdynamics(lattice, model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end
	mpsK = complex(mpsK)
	bath2 = bcsbath(fermionicbath(semicircular(t=1), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
	cache = environments(lattice, mpsK, mpsI)
	g1 = cached_Gτ(lattice, mpsK, mpsI, cache=cache)


	disbath2 = discretebath(bath2, δw=0.01)
	ed_model = Toulouse(disbath2, ϵ_d=-ϵ_d)

	τs = collect(0:δτ:β)
	g2 = toulouse_Gτ(ed_model, τs)

	@test norm(g1-g2) / norm(g2) < rtol
end


@testset "BCS Imaginary time evolution 2" begin
	tol = 1.0e-5
	rtol = 1.0e-2
	δτ=0.1
	Nτ = 10
	β = Nτ * δτ
	chi = 60

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10)

	U = 1
	# ϵ_d = -0.7
	ϵ_d = 0.8
	ω = 1
	α = 0.5

	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, bands=2)
	model = bcs_siam(U=U, μ=-ϵ_d)
	mpsK = accsysdynamics(lattice, model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end

	# compare with ED
	Δ = 0.6
	bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc)
	mpsI′ = hybriddynamics!(vacuumstate(lattice), lattice, corr, orbital=1, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) < tol
	cache = environments(lattice, mpsK, mpsI)
	g1 = cached_Gτ(lattice, mpsK, mpsI, cache=cache)

	H, a, adag, H0 = bcs_operators2(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)

	g2 = correlation_2op_1τ(H, a, adag, 0:δτ:β, β=β)

	@test norm(g1-g2) / norm(g2) < rtol

	# complex gap
	Δ = 0.2 + 0.45*im
	mpsK = complex(mpsK)
	bath2 = bcsbath(fermionicbath(DiracDelta(ω=ω, α=α), β=β), Δ=Δ)
	corr = correlationfunction(bath2, lattice)
	mpsI = hybriddynamics_naive(lattice, corr, orbital=1, trunc=trunc)
	mpsI′ = hybriddynamics(lattice, corr, orbital=1, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) < tol
	cache = environments(lattice, mpsK, mpsI)
	g1 = cached_Gτ(lattice, mpsK, mpsI, cache=cache)

	H, a, adag, H0 = bcs_operators2(U, ϵ_d, ω₀=ω, α=α, Δ=Δ)

	g2 = correlation_2op_1τ(H, a, adag, 0:δτ:β, β=β)

	@test norm(g1-g2) / norm(g2) < rtol

end