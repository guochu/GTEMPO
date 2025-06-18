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


	println("computing MPS-IF...")
	bath = bosonicbath(DiracDelta(ω=ω₀, α=α₀), β=β)
	corr = correlationfunction(bath, flattice)
	@time fmpsI1 = hybriddynamics(flattice, corr, trunc=trunc)

	println("bond dimension of mpsI is ", bond_dimension(fmpsI1))


	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)


	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band, trunc=trunc)
		# mpsI2 = bulkconnection!(mpsI2, lattice, band=band, trunc=trunc)
	end

	# mpsK = systhermalstate!(mpsK, lattice, exact_model, β=β)
	mpsI1 = reweighting!(lattice, mpsK, flattice, fmpsI1, trunc=trunc)

	println("bond dimension of bosonic adt is ", bond_dimension(mpsI1))


	cache = environments(lattice, mpsI1)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time g₁ = [cached_greater(lattice, k, mpsI1, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:N+1]
	@time g₂ = [cached_lesser(lattice, k, mpsI1, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:N+1]

	g₁, g₂ = -im*g₁, -im*g₂

	return g₁, g₂

end