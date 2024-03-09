push!(LOAD_PATH, "../../../src")


using GTEMPO, QuadGK
using JSON, Serialization


function J(D::Real, ω::Real)
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1
end

spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


# function main_analytic(t::Float64; ϵ_d=1., δt=0.05)
# 	D = 2.
# 	N = round(Int, t / δt)
# 	ts = [i*δt for i in 0:N]
# 	g = [G(t, spectrum_func(D), -ϵ_d, 0.) for t in ts]	

# 	g = real(g)
# 	data_path = "result/thouless_analytic_real_t$(t)_mu$(ϵ_d)_dt$(δt).json"

# 	results = Dict("ts"=>ts, "gf" => g)

# 	open(data_path, "w") do f
# 		write(f, JSON.json(results))
# 	end

# 	return g
# end

function main(t::Real; ϵ_d=1., δt=0.05, δτ=0.1, order=8, β=40, chi=1024)
	D = 2.
	β = convert(Float64, β)
	t = convert(Float64, t)
	Nt = round(Int, t / δt)
	Nτ = round(Int, β / δτ)
	println("Nt=", Nt, " δt=", δt, "Nτ=", Nτ, " δτ=", δτ, " ϵ_d=", ϵ_d, " β=", β, " chi=", chi, " order=", order)

	ts = [i*δt for i in 1:Nt+1]
	τs = [i*δτ for i in 1:Nτ+1]

	tol = 10.0^(-order)

	trunc = truncdimcutoff(D=chi, ϵ=tol, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, μ = -ϵ_d, U=0.)

	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1)
	println("number of sites, ", length(lattice))

	mpspath = "data/thouless_tempo_mixed_beta$(β)_dtau$(δτ)_t$(t)_dt$(δt)_order$(order)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		corr = correlationfunction(exact_model.bath, lattice)
		@time mpsI = hybriddynamics(lattice, corr, trunc=trunc, band=1)

		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, mpsI)
	end

	println("mpsI bond dimension is ", bond_dimension(mpsI))
	# println("IF Z is ", integrate(mpsI1, lattice), " ", integrate(mpsI2, lattice))
	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	mpsK = boundarycondition(mpsK, lattice, band=1)

	cache = environments(lattice, mpsK, mpsI)

	band = 1

	@time g₁ = [cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=band, cache=cache) for k in 1:Nt+1]
	@time g₂ = [cached_Gm(lattice, 1, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=band, cache=cache) for k in 1:Nt+1]
	@time g₃ = [cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:τ, b2=:τ, band=band, cache=cache) for k in 1:Nτ+1]


	# ts = [i*δt for i in 1:N]

	data_path = "result/thouless_tempo_mixed_beta$(β)_dtau$(δτ)_t$(t)_mu$(ϵ_d)_dt$(δt)_order$(order)_chi$(chi).json"

	results = Dict("ts"=>ts, "taus"=>τs, "bd"=>bond_dimensions(mpsI), "gt"=>g₁, "lt"=>g₂, "gtau"=>g₃)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return real(g₁)
end

function main_vs_D(t::Real; ϵ_d=1., δt=0.05, β=40)
	for chi in 20:20:200
		main(t; ϵ_d=ϵ_d, δt=δt, β=β, chi=chi, order=10)
	end

end
