push!(LOAD_PATH, "../../../src")

using GTEMPO, JSON, Serialization


J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(t; β=1., U=1., ϵ_d=U/2, δt=0.1, chi=60)
	# β = 5.
	D = 1.
	N = round(Int, t / δt)

	ts = [i*δt for i in 0:N]


	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	lattice = GrassmannLattice(N=N, δt=t/N, bands=2, contour=:real)

	println("number of sites ", length(lattice))

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, U=U, μ=-ϵ_d)
	corr = correlationfunction(bath, lattice)

	mpspath = "data/tempo_beta$(β)_t$(t)_dt$(δt)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		lattice2 = similar(lattice, bands=1)
		@time mpsI = hybriddynamics(lattice2, corr, trunc=trunc)
		Serialization.serialize(mpspath, mpsI)
		println("save MPS-IF to path ", mpspath)
	end

	mpsI1 = fillband(lattice, mpsI, band=1)
	mpsI2 = fillband(lattice, mpsI, band=2)

	# println("Z is ", integrate(mpsI, lattice))
	println("bond dimension of mpsI is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))

	@time mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))

	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	println("bond dimension of mpsK is ", bond_dimension(mpsK))

	cache = environments(lattice, mpsK, mpsI1, mpsI2)
	@time gt = [cached_greater(lattice, k, mpsK, mpsI1, mpsI2, cache=cache) for k in 1:lattice.kt]
	@time lt = [cached_lesser(lattice, k, mpsK, mpsI1, mpsI2, cache=cache) for k in 1:lattice.kt]

	data_path = "result/anderson_tempo_beta$(β)_t$(t)_dt$(δt)_U$(U)_e$(ϵ_d)_chi$(chi).json"

	results = Dict("ts"=>ts, "gt" => gt, "lt"=>lt, "bd"=>bond_dimensions(mpsI1))

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return gt, lt

end
