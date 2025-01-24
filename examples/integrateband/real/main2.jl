push!(LOAD_PATH, "../../../src")

using GTEMPO, JSON, Serialization


J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(t; β=1., U=1., ϵ_d=U/2, δt=0.1, chi=60, chi2=500)
	# β = 5.
	D = 1.
	N = round(Int, t / δt)

	ts = [i*δt for i in 0:N]


	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	lattice = GrassmannLattice(N=N, δt=t/N, bands=2, contour=:real)

	println("number of sites ", length(lattice))
	lattice1 = similar(lattice, bands=1)

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, U=U, μ=-ϵ_d)
	corr = correlationfunction(bath, lattice)

	mpspath = "data/tempo_beta$(β)_t$(t)_dt$(δt)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		
		@time mpsI = hybriddynamics(lattice1, corr, trunc=trunc)
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, mpsI)
	end

	# println("Z is ", integrate(mpsI, lattice))
	println("bond dimension of mpsI is ", bond_dimension(mpsI))

	@time mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
	# println("bond dimension of mpsK is ", bond_dimension(mpsK))
	# println("mpsK scale is ", scaling(mpsK))

	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	println("bond dimension of mpsK is ", bond_dimension(mpsK))

	mpsI1 = fillband(lattice, mpsI, band=1)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10, add_back=0)
	mpsKI1 = mpsK * mpsI1
	mpsKI1′ = integrateband(lattice, mpsKI1, band=1)
	canonicalize!(mpsKI1′, alg = Orthogonalize(trunc=trunc2))

	println("bond dimension of mpsKI1 after compression:", bond_dimension(mpsKI1′))

	algmult = DMRGMult1(trunc=trunc2)
	mps_adt = mult(mpsKI1′, mpsI, algmult)
	println("bond dimension of ADT is ", bond_dimension(mps_adt))

	cache = environments(lattice1, mps_adt)
	@time gt = cached_greater_fast(lattice1, mps_adt, cache=cache) 
	@time lt = cached_lesser_fast(lattice1, mps_adt, cache=cache) 

	data_path = "result/anderson_tempo_beta$(β)_t$(t)_dt$(δt)_U$(U)_e$(ϵ_d)_chi$(chi)_chi2$(chi2).json"

	results = Dict("ts"=>ts, "gt" => gt, "lt"=>lt, "bd"=>bond_dimensions(mpsI1))

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return gt, lt

end
