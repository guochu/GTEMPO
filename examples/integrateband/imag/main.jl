push!(LOAD_PATH, "../../../src")

using GTEMPO, JSON, Serialization


J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(β; U=1., ϵ_d=U/2, δτ=0.1, chi=60)
	# β = 5.
	D = 1.
	N = round(Int, β / δτ)

	τs = [i*δτ for i in 0:N]


	trunc = truncdimcutoff(D=chi, ϵ=1.0e-8, add_back=0)
	lattice = GrassmannLattice(N=N, δτ=β/N, bands=2, contour=:imag)

	println("number of sites ", length(lattice))

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, U=U, μ=-ϵ_d)
	corr = correlationfunction(bath, lattice)

	mpspath = "data/tempo_beta$(β)_N$(N)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		@time mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		@time mpsI2 = hybriddynamics(lattice, corr, trunc=trunc, band=2)
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (mpsI1, mpsI2))
	end

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
	@time g = cached_Gτ(lattice, mpsK, mpsI1, mpsI2, cache=cache)

	data_path = "result/anderson_tempo_beta$(β)_U$(U)_e$(ϵ_d)_N$(N)_chi$(chi).json"

	results = Dict("ts"=>τs, "gf" => g, "bd"=>bond_dimensions(mpsI1))

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g

end
