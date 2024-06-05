# push!(LOAD_PATH, "../../../src")

# using GTEMPO
# using GTEMPO: Δiw_to_Cτ

include("../../../src/includes.jl")

using DelimitedFiles, JSON, Serialization, Interpolations

ρ₀(ϵ, D=1.) = sqrt(1-(ϵ/D)^2) * (D/π)


spectrum_func(D=1.) = SpectrumFunction(ω -> 0.1 * ρ₀(ω, D), lb = -D, ub = D)


function main(; β=10., δτ=0.1, t=2., δt=0.05, U=1., ϵ_d=U/2, chi=100)
	Nτ = round(Int, β / δτ)
	Nt = round(Int, t / δt)
	total_t = 2.
	δt = 0.05
	
	Dh = 0.5
	D = 2 * Dh
	μ = 0.

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	mpath = "beta$(round(Int, β))/"

	lattice = GrassmannLattice(Nτ=Nτ, δτ=δτ, δt=δt, Nt=Nt, bands=2, contour=:mixed)

	# initial guess for Δiw
	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, U=U, μ=-ϵ_d)

	corr1 = correlationfunction(bath, lattice)

	# initial guess for Δiw
	lb = -3.
	ub = 3.
	dw = 1.0e-4
	freqs = collect(frequencies(lb=lb, ub=ub, dw=dw))
	Jw =  [toulouse_Jw(spectrum_func(D), ω) for ω in freqs]

	spec = SpectrumFunction(freqs, Jw)
	bath = fermionicbath(spec, β=β, μ=0)
	corr2 = correlationfunction(bath, lattice)


	return corr1, corr2

end

