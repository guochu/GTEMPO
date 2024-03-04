include("util.jl")

using DelimitedFiles
using Serialization



function main(β, δt=0.05, order=7)
	ϵ_d = -1.
	β = convert(Float64, β)
	U = 0.
	D = 2

	# δt = 0.05
	t = 30.
	N = round(Int, t / δt)
	# N = 600
	# t = N * δt

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, μ = ϵ_d, U=U)

	bands = 1
	lattice_o = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=bands)
	println("number of sites, ", length(lattice_o))


	corr = correlationfunction(exact_model.bath, lattice_o)

	trunc = truncdimcutoff(D=512, ϵ=10.0^(-order), add_back=0)
	truncK = truncdimcutoff(D=512, ϵ=1.0e-10, add_back=0)

	lattice = similar(lattice_o, N=0)
	ns = zeros(Float64, N)
	g₁ = zeros(ComplexF64, N)
	l₁ = zeros(ComplexF64, N)
	bds = zeros(Int, N)

	mpsI = vacuumstate(lattice)
	mpsK = vacuumstate(lattice)

	band = 1

	for k in 2:N+1
		println("the $k-th evolution step...")
		lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
		@assert k == timesteps(mpsI, lattice) == lattice.k

		@time mpsI = hybriddynamicsstepper(mpsI, lattice, corr, trunc=trunc)
		@time mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=truncK)
		mpsK′ = boundarycondition(mpsK, lattice)
		# observables
		cache = environments(lattice, mpsK′, mpsI)
		@time ns[k-1] = cached_occupation(lattice, k-1, mpsK′, mpsI, cache=cache, band=band)
		@time g₁[k-1] = cached_Gt(lattice, k, 1, mpsK′, mpsI, cache=cache, c1=false, c2=true, b1=:+, b2=:+, band=band)
		@time l₁[k-1] = cached_Gt(lattice, 1, k, mpsK′, mpsI, cache=cache, c1=true, c2=false, b1=:-, b2=:+, band=band)
		bds[k-1] = bond_dimension(mpsI)


	end

	ts = [i*δt for i in 1:N]

	println(ns)
	writedlm("result/tempo1_beta$(β)_N$(N)_dt$(δt)_mu$(ϵ_d)_order$(order).dat", [ts real(g₁) imag(g₁) real(l₁) imag(l₁) ns bds])

end

