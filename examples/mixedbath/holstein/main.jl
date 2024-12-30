push!(LOAD_PATH, "../../../src")
using GTEMPO


# include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


function main_real_analytic(ϵ_d; t=1, N=100, ω₀=1, α₀=0.5, order=10)
	δt=t / N
	println("N=", N, " t=", t, " ϵ_d=", ϵ_d, " ω₀=", ω₀, " α₀=", α₀, " order=", order)
	ts = [i*δt for i in 1:N+1]

	spec = semicircular(t=t)

	println(holstein_Gw(spec, -20, ϵ_d=-ϵ_d, ω=ω₀, g=α₀, maxiter=order), " ", holstein_Gw(spec, 20, ϵ_d=-ϵ_d, ω=ω₀, g=α₀, maxiter=order))

	g1 = [holstein_Gt(spec, tj, ϵ_d=-ϵ_d, ω=ω₀, g=α₀, maxiter=order, wmin=-20., wmax=20.) for tj in ts]

	results = Dict("ts"=>ts, "gf" => g1)

	data_path = "result/holstein_analytic_betaInf_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_mu$(ϵ_d)_order$(order).json"

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return ts, g1
end


function main_real(ϵ_d; β=Inf, t=1, N=100, ω₀=1, α₀=0.5, chi = 100)
	# ϵ_d = 0.5
	δt=t / N
	println("N=", N, " t=", t, " ϵ_d=", ϵ_d, " β=", β, " ω₀=", ω₀, " α₀=", α₀, " chi=", chi)

	ts = [i*δt for i in 1:N+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1)
	println("number of sites, ", length(lattice))

	# bath = bosonicbath(spectrum_func(), β=β)
	# corr = correlationfunction(bath, lattice)

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	fbath = fermionicbath(semicircular(t=t), β=β, μ=0)

	mpspath = "data/holstein_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		bath = bosonicbath(DiracDelta(ω=ω₀, α=α₀), β=β)
		corr = correlationfunction(bath, lattice)
		@time mpsI1 = retardedinteractdynamics(lattice, corr, trunc=trunc)

		fcorr = correlationfunction(fbath, lattice)
		@time mpsI2 = hybriddynamics(lattice, fcorr, trunc=trunc)

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (mpsI1, mpsI2))
	end

	println("bond dimension of mpsI is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))

	exact_model = SISB(fbath, U=0., μ=-ϵ_d)
	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	# mpsK = systhermalstate!(mpsK, lattice, exact_model, trunc=trunc, δτ=0.001)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))


	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time g = [cached_greater(lattice, k, mpsK, mpsI1, mpsI2, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:N+1]

	data_path = "result/holstein_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_mu$(ϵ_d)_chi$(chi).json"

	results = Dict("ts"=>ts, "bd1"=>bond_dimensions(mpsI1), "bd2"=>bond_dimensions(mpsI2), "gt"=>g)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end


	return g
end

