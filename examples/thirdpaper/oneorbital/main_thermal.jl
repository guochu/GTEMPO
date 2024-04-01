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

function cached_vacuum(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            cache=environments(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:+, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:+, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    return expectationvalue(t, cache)
end

function cached_nup(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            cache=environments(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:-, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:+, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    return expectationvalue(t, cache)
end

function cached_ndown(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            cache=environments(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:+, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:-, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    return expectationvalue(t, cache)
end

function cached_nn(lattice::RealGrassmannLattice, i::Int, A::GrassmannMPS, B::GrassmannMPS...; 
            cache=environments(lattice, A, B...))
    pos1, pos2 = index(lattice, i, conj=false, branch=:+, band=1), index(lattice, i, conj=true, branch=:-, band=1)
    pos1′, pos2′ = index(lattice, i, conj=false, branch=:+, band=2), index(lattice, i, conj=true, branch=:-, band=2)
    t = GTerm(pos1, pos2, pos1′, pos2′, coeff=1)
    return expectationvalue(t, cache)
end


function main(t::Real, t₀::Real=t/2; U=1., ϵ_d=U/2, δt=0.05, β=40, order=6, chi=1024)
	D = 2.
	β = convert(Float64, β)
	t = convert(Float64, t)
	t₀ = convert(Float64, t₀)
	N = round(Int, t / δt)
	N₀ = round(Int, t₀ / δt)
	println("N=", N, " N₀=", N₀, " δt=", δt, " U=", U, " ϵ_d=", ϵ_d," order=", order, " chi=", chi)


	ts = [i*δt for i in 1:N]

	tol = 10.0^(-order)

	trunc = truncdimcutoff(D=chi, ϵ=tol, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	bath = fermionicbath(spectrum_func(D), β=β, μ=0.)
	exact_model = SISB(bath, μ = -ϵ_d, U = U)

	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=2)
	println("number of sites, ", length(lattice))


	mpspath = "data/anderson_tempo1_beta$(β)_t$(t)_dt$(δt)_order$(order)_chi$(chi).mps"
	if ispath(mpspath)
		println("load MPS-IF from path ", mpspath)
		mpsI1, mpsI2 = Serialization.deserialize(mpspath)
	else
		println("computing MPS-IF...")
		corr = correlationfunction(exact_model.bath, lattice)
		@time mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		# @time mpsI2 = hybriddynamics(lattice, corr, trunc=trunc, band=2)
		@time mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)
		# println("Z is ", integrate(mpsI, lattice))
		println("save MPS-IF to path ", mpspath)
		Serialization.serialize(mpspath, (mpsI1, mpsI2))
	end

	println("mpsI bond dimension is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))
	# println("IF Z is ", integrate(mpsI1, lattice), " ", integrate(mpsI2, lattice))

	@time mpsK = accsysdynamics_fast(lattice, exact_model, trunc=truncK, scaling=1)
	mpsK = systhermalstate!(mpsK, lattice, exact_model, trunc=truncK)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	mpsK = boundarycondition(mpsK, lattice, band=1)
	mpsK = boundarycondition(mpsK, lattice, band=2)
	# println("IF K is ", integrate(mpsK, lattice))
	
	gf_ts = collect(N₀:N)
	g₁ = zeros(ComplexF64, length(gf_ts))
	l₁ = zeros(ComplexF64, length(gf_ts))
	cache = environments(lattice, mpsK, mpsI1, mpsI2)

	@time ns = cached_occupation(lattice, mpsK, mpsI1, mpsI2, cache=cache)
	@time nup_ = [cached_nup(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.N]
	@time ndown_ = [cached_ndown(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.N]
	@time nn_ = [cached_nn(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.N]
	@time vacuum_ = [cached_vacuum(lattice, i, mpsK, mpsI1, mpsI2, cache=cache) for i in 1:lattice.N]


	for i in 1:length(gf_ts)
		g₁[i] = cached_Gt(lattice, N₀+i, N₀, mpsK, mpsI1, mpsI2, cache=cache, c1=false, c2=true, b1=:+, b2=:+)
		l₁[i] = cached_Gt(lattice, N₀, N₀+i, mpsK, mpsI1, mpsI2, cache=cache, c1=true, c2=false, b1=:-, b2=:+)
	end


	data_path = "result/anderson_tempo1_beta$(β)_t$(t)_$(t₀)_U$(U)_mu$(ϵ_d)_dt$(δt)_order$(order)_chi$(chi)_2.json"

	results = Dict("ts"=>ts, "ns" => ns, "bd"=>bond_dimensions(mpsI1), "gt"=>g₁, "lt"=>l₁, "gf_ts"=>gf_ts,
					"nup"=>real(nup_), "ndown"=>real(ndown_), "nn"=>real(nn_), "vacuum"=>real(vacuum_))

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return results
end

function main_all_U(t₀::Real, t::Real = t₀ + 20.; δt=0.05, β=40, order=6, chi=1024)
	for U in 0.1:0.1:1.
		main(t, t₀, U=U, δt=δt, β=β, order=order, chi=chi)
	end
	# for U in [0.2, 0.3,0.4, 0.6,0.7,0.8,0.9]
	# 	main(t, t₀, U=U, δt=δt, β=β, order=order, chi=chi)
	# end
end

function main_all_t(; U=1., δt=0.05, β=40, order=6, chi=1024)
	for t₀ in [5., 10., 15., 20., 40., 80.]
		main(t₀+20, t₀, U=U, δt=δt, β=β, order=order, chi=chi)
	end
end
