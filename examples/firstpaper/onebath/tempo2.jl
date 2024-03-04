include("util.jl")

using DelimitedFiles


function main(β, U)
	ϵ_d = -1.
	β = convert(Float64, β)
	U = convert(Float64, U) 
	D = 2

	δt = 0.05
	N = 600
	t = N * δt

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, μ = ϵ_d, U=U)

	bands = (U ≈ 0) ? 1 : 2
	lattice_o = GrassmannLattice(N=N, δt=δt, contour=:real, order=2, bands=bands)
	println("number of sites, ", length(lattice_o))


	corr = correlationfunction(exact_model.bath, lattice_o)

	trunc = truncdimcutoff(D=256, ϵ=1.0e-7, add_back=0)
	truncK = truncdimcutoff(D=256, ϵ=1.0e-10, add_back=0)
	
	lattice = similar(lattice_o, N=0)
	ns = zeros(Float64, N+1)
	g₁ = zeros(ComplexF64, N+1)
	l₁ = zeros(ComplexF64, N+1)
	bds = zeros(Int, N+1)

	mpsI = vacuumstate(lattice)
	mpsK = vacuumstate(lattice)

	band = 1
	ns[1] = occupation(lattice, mpsK, mpsI, band=band)
	g₁[1] = gf(lattice, 1, 1, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=band)
	l₁[1] = gf(lattice, 1, 1, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=band)
	bds[1] = 1

	for k in 2:N+1
		println("the $k-th evolution step...")
		lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
		@assert k == timesteps(mpsI, lattice) == lattice.k
		@time mpsI2 = hybriddynamicsstepper(mpsI, lattice, corr, finalize=true, trunc=trunc)
		bds[k] = bond_dimension(mpsI2)
		@time mpsI = hybriddynamicsstepper(mpsI, lattice, corr, trunc=trunc)

		@time mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=truncK)
		mpsK′ = mpsK
		for band in 1:lattice.bands
			mpsK′ = boundarycondition(mpsK′, lattice, band=band)
		end

		# observables
		cache = environments(lattice, mpsK′, mpsI2)
		@time ns[k] = cached_occupation(lattice, mpsK′, mpsI2, cache=cache, band=band)
		@time g₁[k] = cached_gf(lattice, k, 1, mpsK′, mpsI2, cache=cache, c1=false, c2=true, b1=:+, b2=:+, band=band)
		@time l₁[k] = cached_gf(lattice, 1, k, mpsK′, mpsI2, cache=cache, c1=true, c2=false, b1=:-, b2=:+, band=band)
	end

	ts = [i*δt for i in 0:N]

	println(ns)
	writedlm("result/tempo2_beta$(β)_U$(U)_N$(N)_dt$(δt)_mu$(ϵ_d).dat", [ts real(g₁) imag(g₁) real(l₁) imag(l₁) ns bds])

end
