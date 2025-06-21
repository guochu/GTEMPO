println("------------------------------------")
println("|        Independent Bosons        |")
println("------------------------------------")

@testset "Independent bosons: imaginary time" begin
	rtol = 1.0e-2
	rtol2 = 1.0e-5
	δτ=0.1
	N = 10
	β = N * δτ
	chi = 100

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1)
	flattice = FockLattice(N=N, δτ=δτ, contour=:imag, order=1)

	for ϵ_d in (-0.5, 0, 0.5)
	# for ϵ_d in (0,)
		for spec in (Leggett(d=3, ωc=1), DiracDelta(ω=1, α=0.5))

			bath = bosonicbath(spec, β=β)
			corr = correlationfunction(bath, flattice)

			mpsI = hybriddynamics(flattice, corr, trunc=trunc)
			mpsI′ = hybriddynamics_naive(flattice, corr, trunc=trunc)
			@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

			exact_model = AndersonIM(U=0., μ=-ϵ_d)
			mpsK = sysdynamics(lattice, exact_model, trunc=truncK)

			adt = reweighting!(lattice, mpsK, flattice, mpsI, trunc=trunc)

			for band in 1:lattice.bands
				adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
			end

			cache = environments(lattice, adt)
			g1 = cached_Gτ(lattice, adt, cache=cache)
			g2 = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, Nτ=N)

			@test norm(g1 - g2) / norm(g1) < rtol	
		end	
	end


	ϵ_d = 0.7
	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=2)
	flattice = FockLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=2)
	for U in (0, 1)
		for spec in (Leggett(d=3, ωc=1), DiracDelta(ω=1, α=0.5))

			bath = bosonicbath(spec, β=β)
			corr = correlationfunction(bath, flattice)

			mpsI = hybriddynamics(flattice, corr, trunc=trunc)
			mpsI′ = hybriddynamics_naive(flattice, corr, trunc=trunc)
			@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

			exact_model = AndersonIM(U=U, μ=-ϵ_d)

			mpsK = sysdynamics(lattice, exact_model, trunc=truncK)

			adt = reweighting!(lattice, mpsK, flattice, mpsI, trunc=trunc)

			for band in 1:lattice.bands
				adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
			end

			cache = environments(lattice, adt)
			g1 = cached_Gτ(lattice, adt, cache=cache)
			g2 = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, Nτ=N, U=U, bands=2)

			@test norm(g1 - g2) / norm(g1) < rtol	
		end	
	end
end

@testset "Independent bosons: real time" begin
	rtol = 5*1.0e-2
	rtol2 = 1.0e-3
	β = 0.1
	δt=0.01
	Nt = 10
	t = Nt * δt
	chi = 60

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1)
	flattice = FockLattice(N=Nt, δt=δt, contour=:real, order=1)
	for ϵ_d in (-0.5, 0, 0.5)

		bath = bosonicbath(DiracDelta(ω=1, α=0.5), β=β)
		corr = correlationfunction(bath, flattice)
		mpsI = hybriddynamics(flattice, corr, trunc=trunc)
		mpsI′ = hybriddynamics_naive(flattice, corr, trunc=trunc)
		@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

		exact_model = AndersonIM(U=0., μ=-ϵ_d)
		mpsK = sysdynamics(lattice, exact_model, trunc=truncK)

		adt = reweighting!(lattice, mpsK, flattice, mpsI, trunc=trunc)
		for band in 1:lattice.bands
			adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
		end

		adt = systhermalstate!(adt, lattice, exact_model, trunc=trunc, β=β)
		cache = environments(lattice, adt)

		g1 = [-im*cached_greater(lattice, k, adt, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
		g2 = [-im*cached_lesser(lattice, k, adt, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]

		# g1′, g2′ = noninteracting_real(ϵ_d)

		H, a, adag, H0 = noninteracting_operators(ϵ_d, ω₀=1, α=0.5, d=100)
		# ρimp = zeros(2,2)
		# ρimp[1,1] = 1
		# ρ = kron(ρimp, exp(-β*H0))
		ρ = exp(-β*H0)
		cache = eigencache(H)
		g1′ = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, cache, reverse = false)
		g2′ = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, cache, reverse = true)

		@test norm(g1 - g1′) / norm(g1) < rtol	
		@test norm(g2 - g2′) / norm(g2) < rtol	
	end


	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1, bands=2)
	flattice = FockLattice(N=Nt, δt=δt, contour=:real, order=1, bands=2)
	ϵ_d = 0.3
	for U in (0, 1)

		bath = bosonicbath(DiracDelta(ω=1, α=0.5), β=β)
		corr = correlationfunction(bath, lattice)
		mpsI = hybriddynamics(flattice, corr, trunc=trunc)
		mpsI′ = hybriddynamics_naive(flattice, corr, trunc=trunc)
		# @test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

		exact_model = AndersonIM(U=U, μ=-ϵ_d)
		mpsK = sysdynamics(lattice, exact_model, trunc=truncK)

		mpsK = systhermalstate!(mpsK, lattice, exact_model, trunc=trunc, β=β)

		for band in 1:lattice.bands
			mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
		end

		adt = reweighting!(lattice, mpsK, flattice, mpsI, trunc=trunc)
		
		cache = environments(lattice, adt)

		g1 = [-im*cached_greater(lattice, k, adt, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
		g2 = [-im*cached_lesser(lattice, k, adt, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]

		# g1′, g2′ = interacting_real(U, ϵ_d)
		H, a, adag, H0 = interacting_operators(U, ϵ_d, ω₀=1, α=0.5, d=100)
		ρ = exp(-β*H0)
		cache = eigencache(H)
		g1′ = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, cache, reverse = false)
		g2′ = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, cache, reverse = true)


		@test norm(g1 - g1′) / norm(g1) < rtol	
		@test norm(g2 - g2′) / norm(g2) < rtol	

	end

end


@testset "Independent bosons: mixed time" begin
	rtol = 5.0e-2
	rtol2 = 1.0e-3
	β = 0.5
	δτ = 0.1
	Nτ = round(Int, β/δτ)
	δt=0.05
	Nt = 5
	t = Nt * δt
	chi = 120

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)


	# noninteracting case
	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1)
	flattice = FockLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1)

	ϵ_d = 0.5
	spec = Leggett(d=3, ωc=1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, flattice)

	mpsI = hybriddynamics(flattice, corr, trunc=trunc)
	mpsI′ = hybriddynamics_naive(flattice, corr, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) <= rtol2

	exact_model = AndersonIM(U=0., μ=-ϵ_d)

	mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	adt = reweighting!(lattice, mpsK, flattice, mpsI, trunc=trunc)
	for band in 1:lattice.bands
		adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
	end

	cache = environments(lattice, adt)

	g1 = [-im*cached_Gm(lattice, k, 1, adt, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	g2 = [im*cached_Gm(lattice, 1, k, adt, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	g3 = [cached_Gm(lattice, k, 1, adt, c1=false, c2=true, b1=:τ, b2=:τ, band=1, cache=cache) for k in 1:Nτ+1]

	g1′ = [independentbosons_greater(spec, tj, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
	g2′ = [independentbosons_lesser(spec, tj, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
	g3′ = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, Nτ=Nτ)

	@test norm(g1 - g1′) / norm(g1) < rtol	
	@test norm(g2 - g2′) / norm(g2) < rtol	
	@test norm(g3 - g3′) / norm(g2) < rtol	

	# interacting case
	flattice = FockLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=2)
	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=2)

	ϵ_d = 0.5
	U = -0.8
	spec = Leggett(d=3, ωc=1)
	bath = bosonicbath(spec, β=β)
	corr = correlationfunction(bath, flattice)

	mpsI = hybriddynamics(flattice, corr, trunc=trunc)
	mpsI′ = hybriddynamics_naive(flattice, corr, trunc=trunc)
	@test distance(mpsI, mpsI′) / norm(mpsI) <= 1.0e-2
	
	exact_model = AndersonIM(U=U, μ=-ϵ_d)

	mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	adt = reweighting!(lattice, mpsK, flattice, mpsI, trunc=trunc)
	for band in 1:lattice.bands
		adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
	end

	cache = environments(lattice, adt)
	
	g1 = [-im*cached_Gm(lattice, k, 1, adt, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	g2 = [im*cached_Gm(lattice, 1, k, adt, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	g3 = [cached_Gm(lattice, k, 1, adt, c1=false, c2=true, b1=:τ, b2=:τ, band=1, cache=cache) for k in 1:Nτ+1]

	g1′ = [independentbosons_greater(spec, tj, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]
	g2′ = [independentbosons_lesser(spec, tj, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]
	g3′ = independentbosons_Gτ(spec, β=β, ϵ_d=-ϵ_d, Nτ=Nτ, U=U, bands=2)

	@test norm(g1 - g1′) / norm(g1) < rtol	
	@test norm(g2 - g2′) / norm(g2) < rtol	
	@test norm(g3 - g3′) / norm(g2) < rtol	

end