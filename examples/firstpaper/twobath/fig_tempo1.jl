### Fig.3 of "Efficient method for quantum impurity problems out of equilibrium"

push!(LOAD_PATH, "../../../src")

using GTEMPO
using DelimitedFiles, JSON, Serialization

function J(D::Real, ω::Real)
	# t′ = 0.3162 
	# return (sqrt(4*t^2-ω^2) / (2*π*t^2)) * t1^2 
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1 / 2
	# return (2/(D*pi)) * sqrt(1 - (ω/D)^2 ) * t′^2 
end


spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)

# Vs = [0.       0.35714286  0.71428571  1.07142857  1.42857143  1.78571429, 2.14285714  2.5         2.85714286]

const Vs = [0., 0.17857143, 0.35714286, 0.53571429, 0.71428571, 0.89285715, 1.07142857, 1.25, 1.42857143, 1.60714287, 1.78571429, 1.96428573, 2.14285714, 2.32142859, 2.5, 2.67857145, 2.85714286]

function main_tempo_1order_b_all(V_over_Gamma, tall=4.2, dt =0.007, order=7, chi=60)
	for U in [0., 2., 4., 6., 8.]
		main_tempo_1order_b(V_over_Gamma, U, tall, dt, order, chi)
	end
end

function main_tempo_1order_b(V_over_Gamma=0.1, U_over_Gamma=0., tall=4.2, dt=0.007, order=6, chi=60)
	println("run for V = ", V_over_Gamma, ", U = ", U_over_Gamma)
	β = Inf
	D = 2

	Γ = 0.1
	V = V_over_Gamma * Γ
	U = U_over_Gamma * Γ

	leftmu = V / 2
	rightmu = -V / 2
	ϵ_d = 0

	δt = dt / Γ
	# t = 4.2 / Γ
	t = tall / Γ
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


	trunc = truncdimcutoff(D=chi, ϵ=10.0^(-order), add_back=0)
	mpspath = "data/tempo1_beta$(β)_N$(N)_V$(V_over_Gamma)_dt$(dt)_order$(order)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		corr = leftcorr + rightcorr
		@time mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		# @time mpsI2 = hybriddynamics(lattice, corr, trunc=trunc, band=2)
		@time mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)
		# println("Z is ", integrate(mpsI, lattice))
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (mpsI1, mpsI2))
	end
	bds = bond_dimensions(mpsI1)

	println("mpsI bond dimension is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))
	# println("IF Z is ", integrate(mpsI1, lattice), " ", integrate(mpsI2, lattice))

	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	@time mpsK = sysdynamics!(vacuumstate(lattice), lattice, exact_model, trunc=truncK)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	mpsK = boundarycondition(mpsK, lattice, band=1)
	mpsK = boundarycondition(mpsK, lattice, band=2)
	# println("IF K is ", integrate(mpsK, lattice))
	
	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	@time ns2 = cached_occupation(lattice, mpsK, mpsI1, mpsI2, cache=cache)

	N_s = round(Int, 4.2 / dt)
	@time currents_left2 = [cached_electriccurrent_fast(lattice, leftcorr, k+1, mpsK, mpsI1, mpsI2, cache=cache) for k in N_s:10:N]
	@time currents_right2 = [cached_electriccurrent_fast(lattice, rightcorr, k+1, mpsK, mpsI1, mpsI2, cache=cache) for k in N_s:10:N]
	currents_ts = collect(N_s:10:N) * dt

	data_path = "result/anderson_tempo1_V$(V_over_Gamma)_U$(U_over_Gamma)_N$(N)_dt$(dt)_order$(order)_chi$(chi).json"

	results = Dict("ts"=>ts, "ns"=>ns2, "Ileft" => real(currents_left2), "ts_I"=>currents_ts, "Iright"=>real(currents_right2), "bd"=>bds)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end
	# return ts, ns2, currents_left2, currents_right2
end
