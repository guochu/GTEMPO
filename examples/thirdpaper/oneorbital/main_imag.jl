push!(LOAD_PATH, "../../../src")

using GTEMPO, QuadGK
using JSON, Serialization


function J(D::Real, ω::Real)
	# t′ = 0.3162 
	# return (sqrt(4*t^2-ω^2) / (2*π*t^2)) * t1^2 
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1 
	# return (2/(D*pi)) * sqrt(1 - (ω/D)^2 ) * t′^2 
end

spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(β; U=1., ϵ_d=U/2, δτ=0.1, order=10, chi=1024)
	# β = 5.
	D = 2.
	N = round(Int, β / δτ)

	τs = [i*δτ for i in 0:N]

	tol = 10.0^(-order)
	trunc = truncdimcutoff(D=chi, ϵ=tol, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	lattice = GrassmannLattice(N=N, δτ=β/N, bands=2, contour=:imag)

	println("number of sites ", length(lattice))

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = AndersonIM(U=U, μ=-ϵ_d)
	corr = correlationfunction(bath, lattice)

	mpspath = "data/andersontempo1_beta$(β)_dt$(δτ)_imag_order$(order)_chi$(chi).mps"
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

	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))
	mpsKs = [mpsK]
	for band in 1:lattice.bands
		mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
	end
	println("bond dimension of mpsK is ", bond_dimension(mpsKs[1]), ", number of Ks ", length(mpsKs))


	@time g = parallel_Gτ(lattice, mpsKs, mpsI1, mpsI2)

	data_path = "result/anderson_tempo1_beta$(β)_U$(U)_mu$(ϵ_d)_N$(N)_dt$(δτ)_imag_order$(order)_chi$(chi).json"

	results = Dict("ts"=>τs, "gf" => g, "bd"=>bond_dimensions(mpsI1))

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g

end

