println("------------------------------------")
println("|              IRLM Model          |")
println("------------------------------------")

function free_irlm(;μ::Real, J::Real)
	hsys = zeros(Float64, 3, 3)
	hsys[2, 2] = μ
	hsys[1, 2] = J
	hsys[2, 1] = J
	hsys[2, 3] = J
	hsys[3, 2] = J
	return hsys
end

@testset "IRLM: real time transport" begin
	rtol = 5.0e-2

	β = 10
	μ = 0.5
	J = 1.
	V = 1.

	N = 20
	δt = 0.02
	t = N*δt
	ts = [i*δt for i in 1:N]

	dw = 0.01

	leftmu = V / 2
	rightmu = -V / 2

	# tol = 1.0e-2

	# ED
	leftbath = fermionicbath(semicircular(), β=β, μ=leftmu)
	rightbath = fermionicbath(semicircular(), β=β, μ=rightmu)

	leftconfig = discretebath(leftbath, δw=dw)
	rightconfig = discretebath(rightbath, δw=dw)

	hsys = free_irlm(μ=μ, J=J)
	model = BoundaryDriving(hsys, leftconfig, rightconfig)

	ham = cmatrix(model)
	observer1 = leftparticlecurrent_cmatrix(model)
	observer2 = rightparticlecurrent_cmatrix(model)

	ρ₀ = separablestate(model, zeros(3,3))

	syssite1, syssite2, syssite3 = 1, 2, 3
	ns1 = Float64[]
	ns2 = Float64[]
	ns3 = Float64[]
	currents1 = ComplexF64[]
	currents2 = ComplexF64[]
	cache = freefermions_cache(ham)
	for i in 1:N
		ρ = freefermions_timeevo(ρ₀, ham, ts[i], cache)
		push!(currents1, sum(observer1 .* ρ))
		push!(currents2, sum(observer2 .* ρ))
		push!(ns1, real(ρ[syssite1, syssite1]))
		push!(ns2, real(ρ[syssite2, syssite2]))
		push!(ns3, real(ρ[syssite3, syssite3]))
	end	
	# currents1 = -2*im .* currents1
	# currents2 = -2*im .* currents2

	# GTEMPO
	lattice = GrassmannLattice(δt=δt, N=N, contour=:real, bands=3)
	leftcorr = correlationfunction(leftbath, lattice)
	rightcorr = correlationfunction(rightbath, lattice)

	chi = 60
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=1000, ϵ=1.0e-10, add_back=0)

	exact_model = IRLM(μ=μ, J=J, U=0)
	mpsK = accsysdynamics_fast(lattice, exact_model, trunc=truncK, scaling=100)
	for band in 1:lattice.bands
		boundarycondition!(mpsK, lattice, trunc=truncK, band=band)
	end

	# corr = leftcorr + rightcorr
	mpsI1 = hybriddynamics(lattice, leftcorr, trunc=trunc, band=1)
	mpsI2 = hybriddynamics(lattice, rightcorr, trunc=trunc, band=3)

	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	ns1′ = cached_occupation(lattice, mpsK, mpsI1, mpsI2, cache=cache, band=1)
	@test norm(ns1 - ns1′) / norm(ns1) < rtol
	ns2′ = cached_occupation(lattice, mpsK, mpsI1, mpsI2, cache=cache, band=2)
	@test norm(ns2 - ns2′) / norm(ns2) < 0.2
	ns3′ = cached_occupation(lattice, mpsK, mpsI1, mpsI2, cache=cache, band=3)
	@test norm(ns3 - ns3′) / norm(ns3) < rtol
	currents1′ = cached_electriccurrent_fast(lattice, leftcorr, mpsK, mpsI1, mpsI2, cache=cache, band=1)
	@test norm(currents1 - currents1′) / norm(currents1) < rtol
	currents2′ = cached_electriccurrent_fast(lattice, rightcorr, mpsK, mpsI1, mpsI2, cache=cache, band=3)
	@test norm(currents2 - currents2′) / norm(currents2) < rtol
end