push!(LOAD_PATH, "../../../src")
using GTEMPO


# include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization

function model(;μ::Real, U::Real, J::Real)
	h = ImpurityHamiltonian(bands=2)
	push!(h, interaction(1,2,2,1, coeff=U))
	for band in 1:h.bands
		push!(h, tunneling(band, band, coeff=μ))
	end
	push!(h, tunneling(1, 2, coeff=J))
	push!(h, tunneling(2, 1, coeff=J))
	return h
end

function main_real(U, J, ϵ_d=U/2; β=1, t=1, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1, chi = 100)
	# ϵ_d = 0.5
	δt=t / N
	println("N=", N, " t=", t, " U=", U, " J=", J, " ϵ_d=", ϵ_d, " β=", β, " ω₀=", ω₀, " α₀=", α₀, " ω₁=", ω₁, " α₁=", α₁, " chi=", chi)

	ts = [i*δt for i in 1:N+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=2)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=N, δt=δt, contour=:real, order=1, bands=2)

	# bath = bosonicbath(spectrum_func(), β=β)
	# corr = correlationfunction(bath, lattice)

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	fbath = fermionicbath(DiracDelta(ω=ω₁, α=α₁), β=β, μ=0)

	# exact_model = model(U=U, μ=-ϵ_d, J=J)
	exact_model = AndersonIM(U=U, μ=-ϵ_d)

	# mpspath = "data/noninteracting_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_omega1$(ω₁)_alpha1$(α₁)_chi$(chi).mps"
	# if ispath(mpspath)
	# 	println("load MPS-IF from path ", mpspath)
	# 	fmpsI1, mpsI2 = Serialization.deserialize(mpspath)
	# else
		println("computing MPS-IF...")
		bath = bosonicbath(DiracDelta(ω=ω₀, α=α₀), β=β)
		corr = correlationfunction(bath, flattice)
		@time fmpsI1 = hybriddynamics(flattice, corr, trunc=trunc)

		fcorr = correlationfunction(fbath, lattice)
		@time mpsI2 = hybriddynamics(lattice, fcorr, band=1, trunc=trunc)
		for band in 1:lattice.bands
			mpsI2 = boundarycondition!(mpsI2, lattice, band=band, trunc=trunc)
			# mpsI2 = bulkconnection!(mpsI2, lattice, band=band, trunc=trunc)
		end

		

		# mpsI2 = systhermalstate!(mpsI2, lattice, exact_model, β=β)

	# 	println("save MPS-IF to path ", mpspath)
	# 	Serialization.serialize(mpspath, (fmpsI1, mpsI2))
	# end

	println("bond dimension of mpsI is ", bond_dimension(fmpsI1), " ", bond_dimension(mpsI2))


	mpsK = accsysdynamics_fast(lattice, exact_model, trunc=trunc)
	mpsK = systhermalstate!(mpsK, lattice, exact_model, β=β, δτ=0.001)
	mpsI1 = reweighting!(lattice, mpsK, flattice, fmpsI1, trunc=trunc)

	println("bond dimension of bosonic adt is ", bond_dimension(mpsI1))


	cache = environments(lattice, mpsI1, mpsI2)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	band = 1
	@time g₁ = [cached_greater(lattice, k, mpsI1, mpsI2, c1=false, c2=true, b1=:+, b2=:+, band=band, cache=cache) for k in 1:N+1]
	@time g₂ = [cached_lesser(lattice, k, mpsI1, mpsI2, c1=true, c2=false, b1=:-, b2=:+, band=band, cache=cache) for k in 1:N+1]
	# @time g₃′ = [cached_nn(lattice, i, 1, mpsI1, mpsI2, cache=cache, b1=:+, b2=:-) for i in 1:N]
	# println("start calculating nn...")
	# g₃ = ComplexF64[]
	# pos2 = index(flattice, 1, branch=:-, band=1)
	# ftmp = apply!(NTerm(pos2, coeff=1), copy(fadt))
	# lattice, mpsItmp = focktograssmann(lattice.ordering, flattice, ftmp, trunc=trunc)
	# v = integrate(lattice, mpsItmp, mpsI2) / Zvalue(cache)
	# push!(g₃, v)
	# for i in 2:N
	# 	pos1 = index(flattice, i, branch=:+, band=1)
	# 	ftmp = apply!(NTerm(pos1, pos2, coeff=1), copy(fadt))
	# 	lattice, mpsItmp = focktograssmann(lattice.ordering, flattice, ftmp, trunc=trunc)
	# 	v = integrate(lattice, mpsItmp, mpsI2) / Zvalue(cache)
	# 	push!(g₃, v)
	# end
	# println("finish calculating nn...")

	# g₁, g₂ = -im*g₁, -im*g₂

	# data_path = "result/noninteracting_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_omega1$(ω₁)_alpha1$(α₁)_mu$(ϵ_d)_chi$(chi).json"

	# results = Dict("ts"=>ts, "bd1"=>bond_dimensions(mpsI1), "bd2"=>bond_dimensions(mpsI2), "gt"=>g₁, "lt"=>g₂, "nn"=>g₃, "nn2"=>g₃′)

	# println("save results to ", data_path)

	# open(data_path, "w") do f
	# 	write(f, JSON.json(results))
	# end


	# return g₁, g₂, g₃

	g₁, g₂ = -im*g₁, -im*g₂
	return g₁, g₂
end