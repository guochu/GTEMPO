push!(LOAD_PATH, "../../../src")
using GTEMPO


# include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


function main_real_analytic(ϵ_d; t=1, N=100, ω₀=1, α₀=0.5, order=10, wmax=100.)
	δt=t / N
	println("N=", N, " t=", t, " ϵ_d=", ϵ_d, " ω₀=", ω₀, " α₀=", α₀, " order=", order)
	ts = [i*δt for i in 1:N+1]

	α = sqrt(α₀)

	spec = semicircular(t=1)

	println(holstein_Gw(spec, -wmax, ϵ_d=-ϵ_d, ω=ω₀, g=α, maxiter=order), " ", holstein_Gw(spec, wmax, ϵ_d=-ϵ_d, ω=ω₀, g=α, maxiter=order))

	g1 = [holstein_Gt(spec, tj, ϵ_d=-ϵ_d, ω=ω₀, g=α, maxiter=order, wmin=-wmax, wmax=wmax) for tj in ts]

	results = Dict("ts"=>ts, "gf" => g1)

	data_path = "result/holstein_analytic_betaInf_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_mu$(ϵ_d)_order$(order)_wmax$(wmax).json"

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return ts, g1
end


function main_imag(ϵ_d; β=5, N=50, ω₀=1, α₀=0.5, chi = 100)
	# ϵ_d = 0.5
	Nτ = N
	δτ = β / Nτ
	println(" ϵ_d=", ϵ_d, " β=", β, " N=", N, " ω₀=", ω₀, " α₀=", α₀, " chi=", chi)

	τs = [i*δτ for i in 1:Nτ+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-14, add_back=0)

	
	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, order=1)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=Nτ, δτ=δτ, contour=:imag, order=1)

	# bath = bosonicbath(spectrum_func(), β=β)
	# corr = correlationfunction(bath, lattice)

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	fbath = fermionicbath(semicircular(t=1), β=β, μ=0)

	mpspath = "data/holstein_imaggtempo_beta$(β)_dtau$(δτ)_omega0$(ω₀)_alpha0$(α₀)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		fmpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		bath = bosonicbath(DiracDelta(ω=ω₀, α=α₀), β=β)
		corr = correlationfunction(bath, flattice)
		@time fmpsI1 = hybriddynamics(flattice, corr, trunc=trunc)

		fcorr = correlationfunction(fbath, lattice)
		@time mpsI2 = hybriddynamics(lattice, fcorr, trunc=trunc)

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (fmpsI1, mpsI2))
	end

	println("bond dimension of mpsI is ", bond_dimension(fmpsI1), " ", bond_dimension(mpsI2))

	exact_model = AndersonIM(U=0., μ=-ϵ_d)

	fadt = sysdynamics!(fmpsI1, flattice, exact_model, trunc=trunc)
	lattice, mpsI1 = focktograssmann(lattice.ordering, flattice, fadt, trunc=trunc)

	for band in 1:lattice.bands
		mpsI1 = boundarycondition!(mpsI1, lattice, band=band, trunc=trunc)
		mpsI1 = bulkconnection!(mpsI1, lattice, band=band, trunc=trunc)
	end
	println("bond dimension of bosonic adt is ", bond_dimension(mpsI1))


	cache = environments(lattice, mpsI1, mpsI2)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time gtau = cached_Gτ_fast(lattice, mpsI1, mpsI2, cache=cache) 
	@time gnn = [cached_nn(lattice, i, 1, mpsI1, mpsI2, cache=cache) for i in 2:N]

	data_path = "result/holstein_imaggtempo_beta$(β)_dtau$(δτ)_omega0$(ω₀)_alpha0$(α₀)_mu$(ϵ_d)_chi$(chi).json"

	results = Dict("taus"=>τs, "bd1"=>bond_dimensions(mpsI1), "bd2"=>bond_dimensions(mpsI2), "gtau"=>gtau, "nn"=>gnn)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end


	return gtau
end


function main_real(ϵ_d; β=Inf, t=1, N=100, ω₀=1, α₀=0.5, chi = 100)
	# ϵ_d = 0.5
	δt=t / N
	println("N=", N, " t=", t, " ϵ_d=", ϵ_d, " β=", β, " ω₀=", ω₀, " α₀=", α₀, " chi=", chi)

	ts = [i*δt for i in 1:N+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)

	
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=N, δt=δt, contour=:real, order=1)

	# bath = bosonicbath(spectrum_func(), β=β)
	# corr = correlationfunction(bath, lattice)

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	fbath = fermionicbath(semicircular(t=1), β=β, μ=0)

	mpspath = "data/holstein_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		fmpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		bath = bosonicbath(DiracDelta(ω=ω₀, α=α₀), β=β)
		corr = correlationfunction(bath, flattice)
		@time fmpsI1 = hybriddynamics(flattice, corr, trunc=trunc)

		fcorr = correlationfunction(fbath, lattice)
		@time mpsI2 = hybriddynamics(lattice, fcorr, trunc=trunc)

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (fmpsI1, mpsI2))
	end

	println("bond dimension of mpsI is ", bond_dimension(fmpsI1), " ", bond_dimension(mpsI2))

	exact_model = AndersonIM(U=0., μ=-ϵ_d)

	fadt = sysdynamics!(fmpsI1, flattice, exact_model, trunc=trunc)
	lattice, mpsI1 = focktograssmann(lattice.ordering, flattice, fadt, trunc=trunc)

	for band in 1:lattice.bands
		mpsI1 = boundarycondition!(mpsI1, lattice, band=band, trunc=trunc)
		mpsI1 = bulkconnection!(mpsI1, lattice, band=band, trunc=trunc)
	end
	println("bond dimension of bosonic adt is ", bond_dimension(mpsI1))


	cache = environments(lattice, mpsI1, mpsI2)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time g = cached_greater_fast(lattice, mpsI1, mpsI2, band=1, cache=cache) 

	data_path = "result/holstein_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_mu$(ϵ_d)_chi$(chi).json"

	results = Dict("ts"=>ts, "bd1"=>bond_dimensions(mpsI1), "bd2"=>bond_dimensions(mpsI2), "gt"=>g)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end


	return g
end

function main_imag_vs_chi(ϵ_d; β=5, N=50, ω₀=1, α₀=0.5)
	for chi in [20, 40,60,80,100, 120, 140]
		main_imag(ϵ_d, β=β, N=N, ω₀=ω₀, α₀=α₀, chi=chi)
	end
end

function main_real_vs_chi_Nt(ϵ_d; t=5, ω₀=1, α₀=0.5)
	for chi in [20, 40,60,80,100]
		for Nt in [50, 100,200,400]
			main_real(ϵ_d, t=t, N=Nt, ω₀=ω₀, α₀=α₀, chi=chi)
		end
	end
end
