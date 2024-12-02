push!(LOAD_PATH, "../../../src")
using GTEMPO


# include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


spectrum_func(;d=3) = Leggett(d=d, ωc=1)

# spectrum_func() = DiracDelta(ω=1, α=0.5)

function main_imag_analytic(ϵ_d; β=1, N=20, d=3)
	g = independentbosons_Gτ(spectrum_func(d=d), β=β, ϵ_d=-ϵ_d, N=N)

	data_path = "result/noninteracting_analytic_imag_beta$(β)_mu$(ϵ_d)_N$(N)_d$(d).json"

	δτ = β / N
	τs = [i*δτ for i in 0:N]

	results = Dict("ts"=>τs, "gf" => g)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g
end

function main_real_analytic(ϵ_d; β=1, t=1, N=100, d=3)
	# ϵ_d = 0.5
	δt=t/N
	ts = collect(0:δt:t)
	g1 = [independentbosons_greater(spectrum_func(d=d), tj, β=β, ϵ_d=-ϵ_d) for tj in ts]
	g2 = [independentbosons_lesser(spectrum_func(d=d), tj, β=β, ϵ_d=-ϵ_d) for tj in ts]

	data_path = "result/noninteracting_analytic_real_beta$(β)_mu$(ϵ_d)_t$(t)_N$(N)_d$(d).json"

	results = Dict("ts"=>ts, "gt" => g1, "lt"=>g2)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g1, g2
end

function main_mixed(ϵ_d; β=1, Nτ=20, t=1, Nt=100, d=3, chi = 100)
	# ϵ_d = 0.5
	δτ = β / Nτ
	δt = t / Nt

	println("Nt=", Nt, " δt=", δt, "Nτ=", Nτ, " δτ=", δτ, " ϵ_d=", ϵ_d, " β=", β, " chi=", chi)

	ts = [i*δt for i in 1:Nt+1]
	τs = [i*δτ for i in 1:Nτ+1]


	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1)
	println("number of sites, ", length(lattice))

	mpspath = "data/noninteracting_mixedgtempo_beta$(β)_dtau$(δτ)_t$(t)_dt$(δt)_d$(d)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		bath = bosonicbath(spectrum_func(d=d), β=β)
		corr = correlationfunction(bath, lattice)
		@time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, mpsI)
	end

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	println("bond dimension of mpsI is ", bond_dimension(mpsI))

	fbath = fermionicbath(semicircular(), β=β, μ=0)
	exact_model = SISB(fbath, U=0., μ=-ϵ_d)
	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	# mpsK = systhermalstate!(mpsK, lattice, exact_model, trunc=trunc, δτ=0.001)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))


	cache = environments(lattice, mpsK, mpsI)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time g₁ = [cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	@time g₂ = [cached_Gm(lattice, 1, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	@time g₃ = [cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:τ, b2=:τ, band=1, cache=cache) for k in 1:Nτ+1]
	@time ns = cached_occupation(lattice, mpsK, mpsI, cache=cache)

	g₁, g₂ = -im*g₁, im*g₂

	data_path = "result/noninteracting_mixedgtempo_beta$(β)_dtau$(δτ)_t$(t)_dt$(δt)_mu$(ϵ_d)_d$(d)_chi$(chi).json"

	results = Dict("ts"=>ts, "taus"=>τs, "bd"=>bond_dimensions(mpsI), "gt"=>g₁, "lt"=>g₂, "gtau"=>g₃, "ns" => ns)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g₁, g₂, real(g₃), ns
end



