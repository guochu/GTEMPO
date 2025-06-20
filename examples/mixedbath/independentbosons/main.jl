push!(LOAD_PATH, "../../../src")
using GTEMPO


# include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


function main_real(ϵ_d; β=1, t=1, N=100, ω₀=1, α₀=0.5, chi = 100)
	# ϵ_d = 0.5
	δt=t / N
	println("N=", N, " t=", t, " ϵ_d=", ϵ_d, " β=", β, " ω₀=", ω₀, " α₀=", α₀, " chi=", chi)

	ts = [i*δt for i in 1:N+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=N, δt=δt, contour=:real, order=1)

	# bath = bosonicbath(spectrum_func(), β=β)
	# corr = correlationfunction(bath, lattice)

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)


	exact_model = AndersonIM(U=0., μ=-ϵ_d)

	mpspath = "data/noninteracting_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega$(ω₀)_alpha$(α₀)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		fmpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		bath = bosonicbath(DiracDelta(ω=ω₀, α=α₀), β=β)
		corr = correlationfunction(bath, flattice)
		@time fmpsI = hybriddynamics(flattice, corr, trunc=trunc)

		println("bond dimension of mpsI is ", bond_dimension(fmpsI))
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, fmpsI)

	end


	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)


	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
		# mpsI2 = bulkconnection!(mpsI2, lattice, band=band, trunc=trunc)
	end

	mpsK = systhermalstate!(mpsK, lattice, exact_model, β=β, δτ=0.001)
	mpsI1 = reweighting!(lattice, mpsK, flattice, fmpsI, trunc=trunc)

	println("bond dimension of bosonic adt is ", bond_dimension(mpsI1))


	cache = environments(lattice, mpsI1)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time g₁ = [cached_greater(lattice, k, mpsI1, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:N+1]
	@time g₂ = [cached_lesser(lattice, k, mpsI1, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:N+1]

	g₁, g₂ = -im*g₁, -im*g₂

	data_path = "result/noninteracting_realgtempo_beta$(β)_t$(t)_dt$(δt)_mu$(ϵ_d)_omega$(ω₀)_alpha$(α₀)_chi$(chi).json"

	results = Dict("ts"=>ts, "bd"=>bond_dimensions(mpsI1), "gt"=>g₁, "lt"=>g₂)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g₁, g₂

end

function main_real_int(U, ϵ_d=U/2; β=1, t=1, N=100, ω₀=1, α₀=0.5, chi = 100)
	# ϵ_d = 0.5
	δt=t / N
	println("N=", N, " t=", t, " U=", U, " ϵ_d=", ϵ_d, " β=", β, " ω₀=", ω₀, " α₀=", α₀, " ω₁=", " chi=", chi)

	ts = [i*δt for i in 1:N+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=2)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=N, δt=δt, contour=:real, order=1, bands=2)

	# bath = bosonicbath(spectrum_func(), β=β)
	# corr = correlationfunction(bath, lattice)

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)


	exact_model = AndersonIM(U=U, μ=-ϵ_d)

	mpspath = "data/interacting_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega$(ω₀)_alpha$(α₀)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		fmpsI = Serialization.deserialize(mpspath)
	else

		println("computing MPS-IF...")
		bath = bosonicbath(DiracDelta(ω=ω₀, α=α₀), β=β)
		corr = correlationfunction(bath, flattice)
		@time fmpsI = hybriddynamics(flattice, corr, trunc=trunc)

		println("bond dimension of mpsI is ", bond_dimension(fmpsI))
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, fmpsI)

	end

	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)

	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
	end

	# mpsK = systhermalstate!(mpsK, lattice, exact_model, β=β, δτ=0.001)
	mpsI1 = reweighting!(lattice, mpsK, flattice, fmpsI, trunc=trunc)



	println("bond dimension of bosonic adt is ", bond_dimension(mpsI1))


	cache = environments(lattice, mpsI1)


	@time g₁ = [cached_greater(lattice, k, mpsI1, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:N+1]
	@time g₂ = [cached_lesser(lattice, k, mpsI1, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:N+1]


	g₁, g₂ = -im*g₁, -im*g₂


	data_path = "result/interacting_realgtempo_beta$(β)_t$(t)_dt$(δt)_U$(U)_mu$(ϵ_d)_omega$(ω₀)_alpha$(α₀)_chi$(chi).json"

	results = Dict("ts"=>ts, "bd"=>bond_dimensions(mpsI1), "gt"=>g₁, "lt"=>g₂)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g₁, g₂

end