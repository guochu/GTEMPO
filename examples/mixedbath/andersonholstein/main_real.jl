push!(LOAD_PATH, "../../../src")
using GTEMPO


# include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization

boson_spectrum_func(;α=1, d=3) = Leggett(d=d, ωc=5, α=α)

function main_real(ϵ_d; β=5, t=5, N=50, d=3, α=1, chi = 100)
	# ϵ_d = 0.5
	δt=t / N
	println(" ϵ_d=", ϵ_d, " β=", β, " t= ", t, " N=", N, " d=", d, " α=", α, " chi=", chi)

	ts = [i*δt for i in 1:N+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-14, add_back=0)

	
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=N, δt=δt, contour=:real, order=1)

	# bath = bosonicbath(spectrum_func(), β=β)
	# corr = correlationfunction(bath, lattice)

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	fbath = fermionicbath(semicircular(t=1), β=β, μ=0)
	bbath = bosonicbath(boson_spectrum_func(d=d, α=α), β=β)

	mpspath = "data/andersonholstein_realgtempo_beta$(β)_t$(t)_dt$(δt)_d$(d)_alpha$(α)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		fmpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		corr = correlationfunction(bbath, flattice)
		@time fmpsI1 = hybriddynamics(flattice, corr, trunc=trunc)

		fcorr = correlationfunction(fbath, lattice)
		@time mpsI2 = hybriddynamics(lattice, fcorr, trunc=trunc)

		for band in 1:lattice.bands
			mpsI2 = boundarycondition!(mpsI2, lattice, band=band, trunc=trunc)
			mpsI2 = bulkconnection!(mpsI2, lattice, band=band, trunc=trunc)
		end

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (fmpsI1, mpsI2))
	end

	println("bond dimension of mpsI is ", bond_dimension(fmpsI1), " ", bond_dimension(mpsI2))

	exact_model = AndersonIM(U=0., μ=-ϵ_d)

	fadt = sysdynamics!(fmpsI1, flattice, exact_model, trunc=trunc)
	lattice, mpsI1 = focktograssmann(lattice.ordering, flattice, fadt, trunc=trunc)

	println("bond dimension of bosonic adt is ", bond_dimension(mpsI1))


	cache = environments(lattice, mpsI1, mpsI2)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time gt = cached_greater_fast(lattice, mpsI1, mpsI2, b1=:+, b2=:+, cache=cache) 
	@time lt = cached_lesser_fast(lattice, mpsI1, mpsI2, b1=:-, b2=:+, cache=cache) 
	@time gnn = [cached_nn(lattice, i, 1, mpsI1, mpsI2, b1=:+, b2=:-, cache=cache) for i in 1:N]

	g₃ = ComplexF64[]
	pos2 = index(flattice, 1, branch=:-, band=1)
	ftmp = apply!(NTerm(pos2, coeff=1), copy(fadt))
	lattice, mpsItmp = focktograssmann(lattice.ordering, flattice, ftmp, trunc=trunc)
	v = integrate(lattice, mpsItmp, mpsI2) / Zvalue(cache)
	push!(g₃, v)
	@time for i in 2:N
		pos1 = index(flattice, i, branch=:+, band=1)
		ftmp = apply!(NTerm(pos1, pos2, coeff=1), copy(fadt))
		lattice, mpsItmp = focktograssmann(lattice.ordering, flattice, ftmp, trunc=trunc)
		v = integrate(lattice, mpsItmp, mpsI2) / Zvalue(cache)
		push!(g₃, v)
	end

	data_path = "result/andersonholstein_realgtempo_beta$(β)_t$(t)_dt$(δt)_d$(d)_alpha$(α)_mu$(ϵ_d)_chi$(chi).json"

	results = Dict("ts"=>ts, "bd1"=>bond_dimensions(mpsI1), "bd2"=>bond_dimensions(mpsI2), "gt"=>gt, "lt"=>lt, "nn"=>g₃, "nn2"=>gnn)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end


end



function main_real_vs_chi(ϵ_d; β=5, t=5, N=50, d=3, α=1)
	for chi in [40,80,120, 160, 200]
		main_real(ϵ_d, β=β, t=t, N=N, d=d, α=α, chi=chi)
	end
end

