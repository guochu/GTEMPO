push!(LOAD_PATH, "../../../src")

using GTEMPO
using DelimitedFiles, JSON, Serialization

function J(D::Real, ω::Real)
	# t′ = 0.3162 
	# return (sqrt(4*t^2-ω^2) / (2*π*t^2)) * t1^2 
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1 
	# return (2/(D*pi)) * sqrt(1 - (ω/D)^2 ) * t′^2 
end

spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(t::Real; U=1., ϵ_d=U/2, δt=0.05, β=40, δτ=0.1, order=10, chi=100)
	D = 2.
	β = convert(Float64, β)
	t = convert(Float64, t)
	Nt = round(Int, t / δt)
	Nτ = round(Int, β / δτ)
	println("Nt=", Nt, " δt=", δt, " Nτ=", Nτ, " δτ=", δτ, " β=", β, " U=", U, " ϵ_d=", ϵ_d," order=", order, " chi=", chi)


	ts = [i*δt for i in 1:Nt+1]
	τs = [i*δτ for i in 1:Nτ+1]

	tol = 10.0^(-order)
	trunc = truncdimcutoff(D=chi, ϵ=tol, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	bath = fermionicbath(spectrum_func(D), β=β, μ=0.)
	exact_model = SISB(bath, μ = -ϵ_d, U = U)

	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1, bands=2)
	println("number of sites, ", length(lattice))

	mpspath = "data/mixed_beta$(β)_dtau$(δτ)_t$(t)_dt$(δt)_order$(order)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		corr = correlationfunction(exact_model.bath, lattice)
		@time mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		@time mpsI2 = hybriddynamics(lattice, corr, trunc=trunc, band=2)
		# println("Z is ", integrate(mpsI, lattice))
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (mpsI1, mpsI2))
	end

	println("mpsI bond dimension is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))
	# println("IF Z is ", integrate(mpsI1, lattice), " ", integrate(mpsI2, lattice))

	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	mpsK = boundarycondition(mpsK, lattice, band=1)
	mpsK = boundarycondition(mpsK, lattice, band=2)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	# println("IF K is ", integrate(mpsK, lattice))
	
	cache = environments(lattice, mpsK, mpsI1, mpsI2)

	band = 1
	@time g₁ = [cached_Gm(lattice, k, 1, mpsK, mpsI1, mpsI2, c1=false, c2=true, b1=:+, b2=:+, band=band, cache=cache) for k in 1:Nt+1]
	@time g₂ = [cached_Gm(lattice, 1, k, mpsK, mpsI1, mpsI2, c1=true, c2=false, b1=:-, b2=:+, band=band, cache=cache) for k in 1:Nt+1]
	@time g₃ = [cached_Gm(lattice, k, 1, mpsK, mpsI1, mpsI2, c1=false, c2=true, b1=:τ, b2=:τ, band=band, cache=cache) for k in 1:Nτ+1]



	data_path = "result/mixed_beta$(β)_dtau$(δτ)_t$(t)_dt$(δt)_U$(U)_mu$(ϵ_d)_order$(order)_chi$(chi).json"

	results = Dict("ts"=>ts, "taus"=>τs, "bd"=>bond_dimensions(mpsI1), "gt"=>g₁, "lt"=>g₂, "gtau"=>g₃)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return results
end

function main_all_U(t::Real; δt=0.05, β=40, δτ=0.1, order=10, chi=100)
	for U in [0.1, 0.5, 1.]
		main(t, U=U, δt=δt, β=β, δτ=δτ, order=order, chi=chi)
	end
end

