push!(LOAD_PATH, "../../../src")

using GTEMPO, QuadGK
using JSON, Serialization


J(D, ε) = sqrt(D^2-ε^2)/pi
spectrum_func(D=1) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)



function Gw(ω::Float64, spectrum, ϵ_d, μ=0)
	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
    1.0/(im*ω-ϵ_d-quadgk(ε -> f(ε)/(im*ω+μ-ε), lb, ub)[1])
end

function G(τ::Float64, spectrum, β, ϵ_d, μ=0)
    res = 0.0
    for n = -100:101
        ω = (2n-1)*π/β
        res += (Gw(ω, spectrum, ϵ_d, μ)-1/(im*ω))*exp(-im*τ*ω)
    end
    res = -(res/β-0.5)
end

function main_analytic(β)
	# β = 1.
	U = 0.
	ϵ_d = -U/2
	D = 1.
	δτ = 1 / 32
	N = round(Int, β / δτ)
	τs = [i*δτ for i in 0:N]
	g = [G(τ, spectrum_func(D), β, ϵ_d, 0.) for τ in τs]	

	g = real(g)
	data_path = "result/anderson_analytic_beta$(β)_N$(N).json"

	results = Dict("ts"=>τs, "gf" => g)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g
end

function main(β; U=1., ϵ_d=U/2, δτ=0.1)
	# β = 5.
	D = 1.
	N = round(Int, β / δτ)

	τs = [i*δτ for i in 0:N]
	chi = 512


	trunc = truncdimcutoff(D=chi, ϵ=1.0e-8, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	lattice = GrassmannLattice(N=N, δτ=β/N, bands=2, contour=:imag)

	println("number of sites ", length(lattice))

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, U=U, μ=-ϵ_d)
	corr = correlationfunction(bath, lattice)

	mpspath = "data/tempo1_beta$(β)_N$(N)_chi$(chi).mps"
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

	data_path = "result/anderson_tempo1_beta$(β)_U$(U)_e$(ϵ_d)_N$(N)_chi$(chi).json"

	results = Dict("ts"=>τs, "gf" => g, "bd"=>bond_dimensions(mpsI1))

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g

end
