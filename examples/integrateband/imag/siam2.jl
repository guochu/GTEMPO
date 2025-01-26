push!(LOAD_PATH, "../../../src")

using GTEMPO, JSON, Serialization

J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(β; U=1., ϵ_d=U/2, δτ=0.1, chi=60, chi2=500)
	# β = 5.
	D = 1.
	N = round(Int, β / δτ)

	τs = [i*δτ for i in 0:N]

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)


	lattice = GrassmannLattice(N=N, δτ=β/N, bands=2, contour=:imag)
	lattice1 = similar(lattice, bands=1)

	println("number of sites ", length(lattice))

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = AndersonIM(U=U, μ=-ϵ_d)
	corr = correlationfunction(bath, lattice)

	mpspath = "data/tempo_beta$(β)_N$(N)_chi$(chi)_1.mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")

		@time mpsI = hybriddynamics(lattice1, corr, trunc=trunc, band=1)
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, mpsI)
	end

	# println("Z is ", integrate(mpsI, lattice))
	println("bond dimension of mpsI is ", bond_dimension(mpsI))

	@time mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	println("bond dimension of mpsK is ", bond_dimension(mpsK))

	mpsI1 = fillband(lattice, mpsI, band=1)
	trunc2 = truncdimcutoff(D=chi2, ϵ=1.0e-10, add_back=0)
	mpsKI1 = mult(mpsK, mpsI1, trunc=trunc2)
	mpsKI1′ = integrateband(lattice, mpsKI1, band=1)
	canonicalize!(mpsKI1′, alg = Orthogonalize(trunc=trunc2))

	mps_adt = mult(mpsKI1′, mpsI, trunc=trunc2)
	println("bond dimension of ADT is ", bond_dimension(mps_adt))

	cache = environments(lattice1, mps_adt)

	@time gtau = cached_Gτ(lattice1, mps_adt, cache=cache)

	data_path = "result/anderson_tempo_beta$(β)_U$(U)_e$(ϵ_d)_N$(N)_chi$(chi)_chi2$(chi2).json"

	results = Dict("ts"=>τs, "gf" => gtau, "bd"=>bond_dimensions(mpsI), "bd2"=>bond_dimension(mps_adt))

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return gtau

end