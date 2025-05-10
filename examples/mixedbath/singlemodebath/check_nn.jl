push!(LOAD_PATH, "../../../src")
using GTEMPO


# include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization
using LinearAlgebra: norm


# spectrum_func() = Leggett(d=3, ωc=1)

# spectrum_func() = DiracDelta(ω=1, α=0.5)

function main_imag(ϵ_d; β=1, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1, chi = 100)
	# ϵ_d = 0.5
	δτ= β / N
	println("N=", N, " ϵ_d=", ϵ_d, " β=", β, " ω₀=", ω₀, " α₀=", α₀, " ω₁=", ω₁, " α₁=", α₁, " chi=", chi)

	τs = [i*δτ for i in 0:N]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, ordering=A1Ā1B1B̄1())
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=N, δτ=δτ, contour=:imag, order=1)

	fbath = fermionicbath(DiracDelta(ω=ω₁, α=α₁), β=β, μ=0)

	mpspath = "data/noninteracting_imaggtempo_beta$(β)_dtau$(δτ)_omega0$(ω₀)_alpha0$(α₀)_omega1$(ω₁)_alpha1$(α₁)_chi$(chi).mps"
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

	@time g₁ = cached_Gτ(lattice, mpsI1, mpsI2, cache=cache)
	@time g₃ = [cached_nn(lattice, i, 1, mpsI1, mpsI2, cache=cache) for i in 2:N]


	data_path = "result/noninteracting_imaggtempo_beta$(β)_dtau$(δτ)_omega0$(ω₀)_alpha0$(α₀)_omega1$(ω₁)_alpha1$(α₁)_mu$(ϵ_d)_chi$(chi).json"

	results = Dict("taus"=>τs, "bd1"=>bond_dimensions(mpsI1), "bd2"=>bond_dimensions(mpsI2), "gf"=>g₁, "nn"=>g₃)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end


	return g₁, g₃
end

function main_real(ϵ_d; β=1, t=1, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1, chi = 100)
	# ϵ_d = 0.5
	δt=t / N
	println("N=", N, " t=", t, " ϵ_d=", ϵ_d, " β=", β, " ω₀=", ω₀, " α₀=", α₀, " ω₁=", ω₁, " α₁=", α₁, " chi=", chi)

	ts = [i*δt for i in 1:N+1]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1)
	println("number of sites, ", length(lattice))
	flattice = FockLattice(N=N, δt=δt, contour=:real, order=1)

	# bath = bosonicbath(spectrum_func(), β=β)
	# corr = correlationfunction(bath, lattice)

	# println("computing MPS-IF...")
	# @time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	fbath = fermionicbath(DiracDelta(ω=ω₁, α=α₁), β=β, μ=0)

	exact_model = AndersonIM(U=0., μ=-ϵ_d)

	mpspath = "data/noninteracting_realgtempo_beta$(β)_t$(t)_dt$(δt)_omega0$(ω₀)_alpha0$(α₀)_omega1$(ω₁)_alpha1$(α₁)_chi$(chi).mps"
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


		for band in 1:lattice.bands
			mpsI2 = boundarycondition!(mpsI2, lattice, band=band, trunc=trunc)
			mpsI2 = bulkconnection!(mpsI2, lattice, band=band, trunc=trunc)
		end
		mpsI2 = systhermalstate!(mpsI2, lattice, exact_model, β=β)

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (fmpsI1, mpsI2))
	end

	println("bond dimension of mpsI is ", bond_dimension(fmpsI1), " ", bond_dimension(mpsI2))


	fadt = sysdynamics!(fmpsI1, flattice, exact_model, trunc=trunc)
	lattice, mpsI1 = focktograssmann(lattice.ordering, flattice, fadt, trunc=trunc)

	println("bond dimension of bosonic adt is ", bond_dimension(mpsI1))


	cache = environments(lattice, mpsI1, mpsI2)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time g₂ = [cached_nn(lattice, i, 1, mpsI1, mpsI2, cache=cache, b1=:+, b2=:-) for i in 2:N]
	println("start calculating nn...")
	g₃ = ComplexF64[]
	pos2 = index(flattice, 1, branch=:-, band=1)
	for i in 2:N
		pos1 = index(flattice, i, branch=:+, band=1)
		ftmp = apply!(NTerm(pos1, pos2, coeff=1), copy(fadt))
		lattice, mpsItmp = focktograssmann(lattice.ordering, flattice, ftmp, trunc=trunc)
		v = integrate(lattice, mpsItmp, mpsI2) / Zvalue(cache)
		push!(g₃, v)
	end
	println("finish calculating nn...")

	println("error between two ways of calculating nn is ", norm(g₃-g₂) / norm(g₂))
end

function main_real_vs_chi(ϵ_d; β=5, t=10, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1)
	for chi in [100, 200, 300, 400, 500, 600, 700]
		main_real(ϵ_d, β=β, t=t, N=N, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁, chi=chi)
	end
end
