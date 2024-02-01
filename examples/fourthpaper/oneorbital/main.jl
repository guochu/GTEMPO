### Fig.3 of "Efficient method for quantum impurity problems out of equilibrium"

push!(LOAD_PATH, "../../../src")

using GTEMPO
using DelimitedFiles, JSON, Serialization

function J(D::Real, ω::Real)
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1 / 2
end


spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)

# Vs = [0.       0.35714286  0.71428571  1.07142857  1.42857143  1.78571429, 2.14285714  2.5         2.85714286]

const Vs = [0., 0.17857143, 0.35714286, 0.53571429, 0.71428571, 0.89285715, 1.07142857, 1.25, 1.42857143, 1.60714287, 1.78571429, 1.96428573, 2.14285714, 2.32142859, 2.5, 2.67857145, 2.85714286]

function main_tempo_1order_b_all(V_over_Gamma, t=6., dt =0.007, order=7, k=5)
	for U in [0., 2., 4., 6., 8.]
		main_tempo_1order_b(V_over_Gamma, U, t, dt, order, k)
	end
end

function main_tempo_1order_b(V_over_Gamma, U_over_Gamma, tall, dt, order, k)
	println("run for V= ", V_over_Gamma, ", U= ", U_over_Gamma, " t= ", tall, " dt= ", dt, " order= ", order, " k= ", k)
	β = Inf
	D = 2
	chi = 160
	prony=1.0e-4

	Γ = 0.1
	V = V_over_Gamma * Γ
	U = U_over_Gamma * Γ

	leftmu = V / 2
	rightmu = -V / 2
	ϵ_d = 0

	δt = dt / Γ
	t = tall / Γ
	# t = 0.07 / Γ
	N = round(Int, t / δt)
	println("total number of steps ", N)
	ts = [i*δt for i in 1:N]
	
	leftbath = fermionicbath(spectrum_func(D), β=β, μ=leftmu)
	rightbath = fermionicbath(spectrum_func(D), β=β, μ=rightmu)
	
	exact_model = SIDB(leftbath, rightbath, μ = ϵ_d - U / 2, U = U)
	bands = 2
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=bands)
	println("number of sites, ", length(lattice))
	leftcorr = correlationfunction(exact_model.leftbath, lattice)
	rightcorr = correlationfunction(exact_model.rightbath, lattice)


	trunc = truncdimcutoff(D=chi, ϵ=10.0^(-order))
	mpspath = "data/ti_N$(N)_V$(V_over_Gamma)_dt$(dt)_order$(order)_chi$(chi)_k$(k).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		corr = leftcorr + rightcorr
		algevo = WII()
		algexpan = PronyExpansion(n=15, tol=prony, verbosity=4)
		algmult = DMRG1(trunc)
		alg = TranslationInvariantIF(algevo=algevo, algexpan=algexpan, algmult=algmult, k=k)

		@time mpsI1 = hybriddynamics(lattice, corr, alg, band=1)
		@time mpsI2 = hybriddynamics(lattice, corr, alg, band=2)
		# println("Z is ", integrate(mpsI, lattice))
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (mpsI1, mpsI2))
	end
	bds = bond_dimensions(mpsI1)

	println("mpsI bond dimension is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))
	# println("IF Z is ", integrate(mpsI1, lattice), " ", integrate(mpsI2, lattice))

	truncK = truncdimcutoff(D=512, ϵ=1.0e-10, add_back=0)
	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	mpsK = boundarycondition!(mpsK, lattice, band=1)
	mpsK = boundarycondition!(mpsK, lattice, band=2)
	# println("IF K is ", integrate(mpsK, lattice))

	println("compute partition function...")
	@time Z = integrate(lattice, mpsK, mpsI1, mpsI2)

	@time currents_left2 = [electriccurrent_fast(lattice, leftcorr, k+1, mpsK, mpsI1, mpsI2, Z=Z) for k in 300:10:N]
	@time currents_right2 = [electriccurrent_fast(lattice, rightcorr, k+1, mpsK, mpsI1, mpsI2, Z=Z) for k in 300:10:N]
	currents_ts = ts[1:10:N]

	data_path = "result/ti_V$(V_over_Gamma)_U$(U_over_Gamma)_N$(N)_dt$(dt)_order$(order)_chi$(chi)_k$(k).json"

	results = Dict("ts"=>ts, "Ileft" => real(currents_left2), "ts_I"=>currents_ts, "Iright"=>real(currents_right2), "bd"=>bds)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end
	# return ts, ns2, currents_left2, currents_right2
end
