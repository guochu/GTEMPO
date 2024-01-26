# push!(LOAD_PATH, "../../../src")
# using GTEMPO

include("../../../src/includes.jl")

using JSON, QuadGK


function J(D::Real, ω::Real)
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1
end

spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)

function Gw(ω::Float64, spectrum, ϵ_d, μ=0)
	δ = 1e-8
	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
	1.0/(ω-ϵ_d-quadgk(ε -> f(ε)/(ω+μ-ε+im*δ), lb, ub)[1])
end

function G(t::Float64, spectrum, ϵ_d, μ=0)
    δ = 1e-8
    A = quadgk(ω -> (Gw(ω, spectrum, ϵ_d, μ)-1.0/(ω+im*δ))*exp(-im*ω*t), -20., 20.)[1]
    im*(A/(2π)-im)
end

function main_analytic(t::Float64; ϵ_d=0, δt=0.05)
	D = 2.
	N = round(Int, t / δt)
	ts = [i*δt for i in 0:N]
	g = [G(t, spectrum_func(D), -ϵ_d, 0.) for t in ts]	

	g = real(g)
	data_path = "result/analytic_mu$(ϵ_d)_dt$(δt)_N$(N).json"

	results = Dict("ts"=>ts, "gf" => g)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g
end

# Γ = 0.1

function main_partial(t; ϵ_d=0, β=20., order=7, δt = 0.05)
	D = 2.
	χ = 50

	N = round(Int, t / δt)

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, μ = -ϵ_d, U=0)

	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=1)
	println("number of sites, ", length(lattice))

	corr = correlationfunction(exact_model.bath, lattice)

	trunc = truncdimcutoff(D=χ, ϵ=10.0^(-order), add_back=0)
	truncK = truncdimcutoff(D=χ, ϵ=1.0e-10, add_back=0)


	println("computing MPS-IF...")
	_t = @elapsed mpsI = hybriddynamics(lattice, corr, trunc=trunc, band=1)

	println("time for building Partial IF ", _t)

	bds = bond_dimensions(mpsI)

	println("mpsI bond dimension is ", bond_dimension(mpsI))


	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	mpsK = boundarycondition(mpsK, lattice, band=1)

	cache = environments(lattice, mpsK, mpsI)

	greater = [cached_gf(lattice, k, 1, mpsK, mpsI, cache=cache, c1=false, c2=true, f1=true, f2=true, band=1) for k in 2:lattice.k]
	lesser = [cached_gf(lattice, 1, k, mpsK, mpsI, cache=cache, c1=true, c2=false, f1=false, f2=true, band=1) for k in 2:lattice.k]
	ns = cached_occupation(lattice, mpsK, mpsI, cache=cache)

	ts = [i*δt for i in 1:N]

	data_path = "result/partial_if_mu$(ϵ_d)_beta$(β)_dt$(δt)_N$(N)_order$(order).json"
	results = Dict("ts"=>ts, "gt"=>greater, "lt"=>lesser, "bd"=>bds, "ns"=>ns, "time"=>_t)
	println("save results to path ", data_path)
	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

end

function main_partial_vs_t()
	main_partial(0.5)
	for t in [10., 20., 30., 40., 50., 60.]
		main_partial(t)
	end
end

function main_partial_vs_order()
	main_partial(0.5)
	t = 30.
	for order in [6,7,8]
		main_partial(t, order=order)
	end
end

function main_ti(t; ϵ_d=0, β = 20., order=7, prony=4, k=5, δt = 0.05)
	D = 2.
	χ = 50

	N = round(Int, t / δt)

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, μ = -ϵ_d, U=0)

	lattice = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=1)
	println("number of sites, ", length(lattice))

	corr = correlationfunction(exact_model.bath, lattice)

	trunc = truncdimcutoff(D=χ, ϵ=10.0^(-order), add_back=0)
	truncK = truncdimcutoff(D=χ, ϵ=1.0e-10, add_back=0)


	println("computing MPS-IF...")
	algevo = WII()
	algexpan = PronyExpansion(n=20, tol=10.0^(-prony), verbosity=4)
	algmult = DMRG1(trunc)
	alg = TranslationInvariantIF(algevo=algevo, algexpan=algexpan, algmult=algmult, k=k)
	_t = @elapsed mpsI = hybriddynamics(lattice, corr, alg, band=1)

	println("time for building Translation Invariant IF ", _t)

	bds = bond_dimensions(mpsI)

	println("mpsI bond dimension is ", bond_dimension(mpsI))


	@time mpsK = accsysdynamics_fast(lattice, exact_model, trunc=truncK)
	println("mpsK bond dimension is ", bond_dimension(mpsK))
	mpsK = boundarycondition(mpsK, lattice, band=1)

	cache = environments(lattice, mpsK, mpsI)

	greater = [cached_gf(lattice, k, 1, mpsK, mpsI, cache=cache, c1=false, c2=true, f1=true, f2=true, band=1) for k in 2:lattice.k]
	lesser = [cached_gf(lattice, 1, k, mpsK, mpsI, cache=cache, c1=true, c2=false, f1=false, f2=true, band=1) for k in 2:lattice.k]
	ns = cached_occupation(lattice, mpsK, mpsI, cache=cache)

	ts = [i*δt for i in 1:N]

	data_path = "result/ti_if_mu$(ϵ_d)_beta$(β)_dt$(δt)_N$(N)_order$(order)_prony$(prony)_k$(k).json"
	println("save results to path ", data_path)
	results = Dict("ts"=>ts, "gt"=>greater, "lt"=>lesser, "bd"=>bds, "ns"=>ns, "time"=>_t)
	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

end


function main_ti_vs_t()
	main_ti(0.5)
	for t in [10., 20., 30., 40., 50., 60.]
		main_ti(t)
	end
end

function main_ti_vs_order()
	main_ti(0.5)
	t = 30.
	for order in [6,7,8]
		main_ti(t, order=order)
	end
end

function main_ti_vs_prony()
	main_ti(0.5)
	t = 30.
	for prony in [3,4]
		main_ti(t, prony=prony)
	end
end

function main_ti_vs_k()
	main_ti(0.5)
	t = 30.
	for k in [5, 10, 15]
		main_ti(t, k=k)
	end
end
