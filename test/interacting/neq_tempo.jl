# TCMPS seems the only way to benchmark tempo for interacting system, but it is too time costy.
println("----------------------------------------------")
println("|     NEQ TEMPO For Interacting System       |")
println("----------------------------------------------")



@testset "TEMPO vs TCMPS: NEQ for single bath" begin

	β = 1
	μ = 0
	ϵ_d = -1
	D = 1
	N = 10
	δt = 0.05
	t = N*δt
	U = 1.
	dw = 0.2
	tol = 1.0e-2

	# TCMPS
	bath = fermionicbath(spectrum_func(D), β=β, μ=μ)
	config = thermofield( chainmapping(star(bath, dw=dw)) )
	model = SISBD(config, U=U, μ=ϵ_d)

	state = complex(separable_state(model, sys_states=[0]))
	trunc = truncdimcutoff(D=100, ϵ=1.0e-6)
	gf_stepper = TEBDStepper(tspan=(0 , -im*δt), stepsize=δt, order=2, trunc=trunc)
	ts = [i * δt for i in 0:N]
	g_tcmps, l_tcmps = gf_greater_lesser_t(model, ts, 1, 1, gf_stepper=gf_stepper, th_state=state)

	@test norm(l_tcmps) < tol

	# tempo
	exact_model = SISB(bath, U=U, μ=ϵ_d)
	lattice = GrassmannLattice(N=N, δt=δt, bands=2, contour=:real, order=1)
	mpsI = hybriddynamics(lattice, exact_model, trunc=trunc) 
	mpsI = boundarycondition(mpsI, lattice, band=1)
	mpsI = boundarycondition(mpsI, lattice, band=2)
	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

	Z = integrate(lattice, mpsK, mpsI)
	g₁ = [-im*gf(lattice, i, 1, mpsK, mpsI, c1=false, c2=true, f1=true, f2=true, band=1, Z=Z) for i in 1:lattice.k]
	l₁ = [-im*gf(lattice, 1, i, mpsK, mpsI, c1=true, c2=false, f1=false, f2=true, band=1, Z=Z) for i in 1:lattice.k]

	g₂ = [-im*gf(lattice, i, 1, mpsK, mpsI, c1=false, c2=true, f1=true, f2=true, band=2, Z=Z) for i in 1:lattice.k]
	l₂ = [-im*gf(lattice, 1, i, mpsK, mpsI, c1=true, c2=false, f1=false, f2=true, band=2, Z=Z) for i in 1:lattice.k]


	@test norm(g_tcmps - g₁) / norm(g_tcmps) < tol
	@test norm(g₁ - g₂) / norm(g₁) < tol
	@test norm(l₁ - l₂) / norm(g₁) < tol

end


@testset "TEMPO vs TCMPS: NEQ for two bath" begin

	epsilon_d = -1.25  * pi
	leftmu = -5
	rightmu = 3
	leftbeta = 1
	rightbeta = 0.5
	U = 3.5*pi
	# U = 0
	Dl = 10.
	Dr = 10.
	rtol = 5*1.0e-2

	δt = 0.02
	N = 5
	t = N*δt
	ts = [i*δt for i in 0:N] 

	leftbath = fermionicbath(spectrum_func(Dl), β=leftbeta, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(Dr), β=rightbeta, μ=rightmu)


	leftdw = 0.4
	rightdw = 0.4
	leftconfig = thermofield(chainmapping(star(leftbath, dw=leftdw)))
	rightconfig = thermofield(chainmapping(star(rightbath, dw=rightdw)))

	# TCMPS
	model = SIDBD(leftconfig, rightconfig, μ=epsilon_d, U=U)
	state = complex(separable_state(model, sys_states=[0]))

	symmetry = SpinCharge()
	sys_site = only(default_sys_sites(model))
	ham = consolidate(FermionicHamiltonian(model), symmetry=symmetry)
	
	observer_n = consolidate(TwoBodyTerm(sys_site, sys_site), symmetry=symmetry) 
	observer1 = consolidate(left_sysbath_tunneling(model), symmetry=symmetry)
	observer2 = consolidate(right_sysbath_tunneling(model), symmetry=symmetry)

	trunc = truncdimcutoff(D=100, ϵ=1.0e-6, add_back=0)

	ns = Float64[]
	currents_left = ComplexF64[]
	currents_right = ComplexF64[]

	local cache
	for i in 1:N
		stepper = TEBDStepper(tspan=(-im*(i-1)*δt, -im*i*δt), stepsize=δt, order=2, trunc=trunc)
		if !(@isdefined cache)
			state, cache = timeevo!(state, ham, stepper)
		else
			state, cache = timeevo!(state, ham, stepper, cache)
		end
		push!(ns, real(expectation(observer_n, state)) )
		push!(currents_left, expectation(observer1, state) )
		push!(currents_right, expectation(observer2, state) )
	end
	# why the convention is so strength for TCMPS?
	currents_left = conj(im .* currents_left)
	currents_right = conj(im .* currents_right)
	ns ./= 2

	# TEMPO
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, bands=2, order=1)
	exact_model = SIDB(leftbath, rightbath, μ=epsilon_d, U=U)
	leftcorr = correlationfunction(exact_model.leftbath, lattice)
	rightcorr = correlationfunction(exact_model.rightbath, lattice)
	mpsI = qim_hybriddynamics(lattice, leftcorr + rightcorr, trunc=trunc) 
	mpsI = boundarycondition(mpsI, lattice, band=1)
	mpsI = boundarycondition(mpsI, lattice, band=2)	
	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
	ns2 = occupation(lattice, mpsK, mpsI)

	currents_left2 = electriccurrent(lattice, leftcorr, mpsK, mpsI)
	currents_right2 = electriccurrent(lattice, rightcorr, mpsK, mpsI)

	@test norm(ns - ns2) / norm(ns) < rtol
	@test norm(currents_left - currents_left2) / norm(currents_left) < rtol
	@test norm(currents_right - currents_right2) / norm(currents_right) < rtol

end