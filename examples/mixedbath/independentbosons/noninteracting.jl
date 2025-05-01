push!(LOAD_PATH, "../../../src")
using GTEMPO


# include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


spectrum_func(;α=1, d=3) = Leggett(d=d, ωc=5, α=α)

# spectrum_func() = DiracDelta(ω=1, α=0.5)

function main_imag_analytic(ϵ_d; β=1, N=20, d=1, α=1)
	g = independentbosons_Gτ(spectrum_func(d=d, α=α), β=β, ϵ_d=-ϵ_d, Nτ=N)

	data_path = "result/noninteracting_analytic_imag_beta$(β)_mu$(ϵ_d)_N$(N)_d$(d)_alpha$(α).json"

	δτ = β / N
	τs = [i*δτ for i in 0:N]

	results = Dict("ts"=>τs, "gf" => g)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g
end

function main_real_analytic(ϵ_d; β=1, t=1, N=100, d=1, α=1)
	# ϵ_d = 0.5
	δt=t/N
	ts = collect(0:δt:t)
	g1 = [independentbosons_greater(spectrum_func(d=d, α=α), tj, β=β, ϵ_d=-ϵ_d) for tj in ts]
	g2 = [independentbosons_lesser(spectrum_func(d=d, α=α), tj, β=β, ϵ_d=-ϵ_d) for tj in ts]

	data_path = "result/noninteracting_analytic_real_beta$(β)_mu$(ϵ_d)_t$(t)_N$(N)_d$(d)_alpha$(α).json"

	results = Dict("ts"=>ts, "gt" => g1, "lt"=>g2)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g1, g2
end

function main_imag(ϵ_d; β=1, Nτ=20, d=1, chi = 100, α=1)
	# ϵ_d = 0.5
	δτ = β / Nτ

	println(" Nτ=", Nτ, " δτ=", δτ, " ϵ_d=", ϵ_d, " β=", β, " chi=", chi, " d=", d, " α=", α)

	τs = [i*δτ for i in 1:Nτ+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(N=Nτ, δτ=δτ, contour=:imag, order=1)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=Nτ, δτ=δτ, contour=:imag, order=1)

	mpspath = "data/noninteracting_imaggtempo_beta$(β)_dtau$(δτ)_d$(d)_alpha$(α)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		fmpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		bath = bosonicbath(spectrum_func(d=d, α=α), β=β)
		corr = correlationfunction(bath, flattice)
		@time fmpsI = hybriddynamics(flattice, corr, trunc=trunc)

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, fmpsI)
	end

	println("bond dimension of fmpsI is ", bond_dimension(fmpsI))

	exact_model = AndersonIM(U=0., μ=-ϵ_d)

	fadt = sysdynamics!(fmpsI, flattice, exact_model, trunc=trunc)
	lattice, adt = focktograssmann(lattice.ordering, flattice, fadt, trunc=trunc)

	
	for band in 1:lattice.bands
		adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
		adt = bulkconnection!(adt, lattice, band=band, trunc=trunc)
	end
	println("bond dimension of bosonic adt is ", bond_dimension(adt))

	cache = environments(lattice, adt)

	@time gtau = [cached_Gτ(lattice, k, 1, adt, c1=false, c2=true, band=1, cache=cache) for k in 1:Nτ+1]

	data_path = "result/noninteracting_imaggtempo_beta$(β)_dtau$(δτ)_mu$(ϵ_d)_d$(d)_alpha$(α)_chi$(chi).json"

	results = Dict("taus"=>τs, "bd"=>bond_dimensions(adt), "gtau"=>gtau)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return τs, gtau
end

function main_mixed(ϵ_d; β=1, Nτ=20, t=1, Nt=100, d=1, chi = 100, α=1)
	# ϵ_d = 0.5
	δτ = β / Nτ
	δt = t / Nt

	println("Nt=", Nt, " δt=", δt, " Nτ=", Nτ, " δτ=", δτ, " ϵ_d=", ϵ_d, " β=", β, " chi=", chi, " d=", d, " α=", α)

	ts = [i*δt for i in 1:Nt+1]
	τs = [i*δτ for i in 1:Nτ+1]


	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1)

	mpspath = "data/noninteracting_mixedgtempo_beta$(β)_dtau$(δτ)_t$(t)_dt$(δt)_d$(d)_alpha$(α)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		fmpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		bath = bosonicbath(spectrum_func(d=d, α=α), β=β)
		corr = correlationfunction(bath, flattice)
		@time fmpsI = hybriddynamics(flattice, corr, trunc=trunc)

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, fmpsI)
	end

	println("bond dimension of fmpsI is ", bond_dimension(fmpsI))

	exact_model = AndersonIM(U=0., μ=-ϵ_d)

	fadt = sysdynamics!(fmpsI, flattice, exact_model, trunc=trunc)
	lattice, adt = focktograssmann(lattice.ordering, flattice, fadt, trunc=trunc)

	
	for band in 1:lattice.bands
		adt = boundarycondition!(adt, lattice, band=band, trunc=trunc)
		adt = bulkconnection!(adt, lattice, band=band, trunc=trunc)
	end
	println("bond dimension of bosonic adt is ", bond_dimension(adt))

	cache = environments(lattice, adt)


	@time g₁ = cached_greater_fast(lattice, adt, band=1, cache=cache) 
	@time g₂ = cached_lesser_fast(lattice, adt, band=1, cache=cache) 
	@time g₃ = cached_Gτ_fast(lattice, adt, band=1, cache=cache) 
	@time ns = cached_occupation(lattice, adt, cache=cache)

	g₁, g₂ = -im*g₁, im*g₂

	data_path = "result/noninteracting_mixedgtempo_beta$(β)_dtau$(δτ)_t$(t)_dt$(δt)_mu$(ϵ_d)_d$(d)_alpha$(α)_chi$(chi).json"

	results = Dict("ts"=>ts, "taus"=>τs, "bd"=>bond_dimensions(adt), "gt"=>g₁, "lt"=>g₂, "gtau"=>g₃, "ns" => ns)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g₁, g₂, real(g₃), ns
end


function main_imag_vs_chi_Ntau(ϵ_d; β=5, d=1, α=1)
	for chi in [20, 40,60,80,100]
		for Ntau in [50, 100,200,400,800]
			main_imag(ϵ_d, β=β, Nτ=Ntau, d=d, α=α, chi=chi)
		end
	end
end

function main_imag_vs_chi(ϵ_d; β=5, Nτ=50, d=1, α=1)
	for chi in [20, 40,60,80,100]
		main_imag(ϵ_d, β=β, Nτ=Nτ, d=d, α=α, chi=chi)
	end
end

function main_mixed_vs_chi_Nt(ϵ_d; β=5, Nτ=50, t=5, d=1, α=1)
	for chi in [60,80,100,120,140]
		for Nt in [50, 100,200,400]
			main_mixed(ϵ_d, β=β, Nτ=Nτ, t=t, Nt=Nt, d=d, α=α, chi=chi)
		end
	end
end


