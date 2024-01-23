push!(LOAD_PATH, "../../../src")


using GTEMPO, QuadGK
using JSON, Serialization


function J(D::Real, ω::Real)
	return (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1
end

spectrum_func(D) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)



# function Gw(ω::Float64, spectrum, ϵ_d, μ=0)
# 	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
#     1.0/(im*ω-ϵ_d-quadgk(ε -> f(ε)/(im*ω+μ-ε), lb, ub)[1])
# end

# function G(τ::Float64, spectrum, β, ϵ_d, μ=0)
#     res = 0.0
#     for n = -100:101
#         ω = (2n-1)*π/β
#         res += (Gw(ω, spectrum, ϵ_d, μ)-1/(im*ω))*exp(-im*τ*ω)
#     end
#     res = -(res/β-0.5)
# end

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

function main_analytic(t::Float64; ϵ_d=1., δt=0.05)
	D = 2.
	N = round(Int, t / δt)
	ts = [i*δt for i in 0:N]
	g = [G(t, spectrum_func(D), -ϵ_d, 0.) for t in ts]	

	g = real(g)
	data_path = "result/thouless_analytic_real_t$(t)_mu$(ϵ_d)_dt$(δt).json"

	results = Dict("ts"=>ts, "gf" => g)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return g
end

function main(t::Real; ϵ_d=1., δt=0.05, β=40)
	D = 2.
	β = convert(Float64, β)
	t = convert(Float64, t)
	N = round(Int, t / δt)
	println("N=", N, " δt=", δt, " ϵ_d=", ϵ_d)

	ts = [i*δt for i in 1:N]
	chi = 1024

	trunc = truncdimcutoff(D=chi, ϵ=1e-6, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, μ = -ϵ_d, U=0.)

	lattice_o = GrassmannLattice(N=N, δt=δt, contour=:real, order=1, bands=1)
	println("number of sites, ", length(lattice_o))

	corr = correlationfunction(exact_model.bath, lattice_o)

	lattice = similar(lattice_o, N=0)
	ns = zeros(Float64, N)
	g₁ = zeros(ComplexF64, N)
	l₁ = zeros(ComplexF64, N)
	bds = zeros(Int, N)

	mpsI = vacuumstate(lattice)
	mpsK = vacuumstate(lattice)

	band = 1

	for k in 2:N+1
		println("the $k-th evolution step...")
		lattice, mpsI, mpsK = makestep(lattice, mpsI, mpsK)
		@assert k == timesteps(mpsI, lattice) == lattice.k

		@time mpsI = hybriddynamicsstepper(mpsI, lattice, corr, trunc=trunc)
		@time mpsK = sysdynamicsstepper!(mpsK, lattice, exact_model, trunc=truncK)
		mpsK′ = boundarycondition(mpsK, lattice, band=band)
		# observables
		cache = environments(lattice, mpsK′, mpsI)
		@time ns[k-1] = cached_occupation(lattice, k-1, mpsK′, mpsI, cache=cache, band=band)
		@time g₁[k-1] = cached_gf(lattice, k, 1, mpsK′, mpsI, cache=cache, c1=false, c2=true, f1=true, f2=true, band=band)
		@time l₁[k-1] = cached_gf(lattice, 1, k, mpsK′, mpsI, cache=cache, c1=true, c2=false, f1=false, f2=true, band=band)
		bds[k-1] = bond_dimension(mpsI)

	end

	# ts = [i*δt for i in 1:N]

	data_path = "result/thouless_tempo_real_beta$(β)_t$(t)_mu$(ϵ_d)_dt$(δt)_e6.json"

	results = Dict("ts"=>ts, "ns" => ns, "bd"=>bds, "gt"=>g₁, "lt"=>l₁)

	println("save results to ", data_path)

	open(data_path, "w") do f
		write(f, JSON.json(results))
	end

	return real(g₁)
end

