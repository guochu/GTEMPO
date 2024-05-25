# push!(LOAD_PATH, "../../../src")

# using GTEMPO
# using GTEMPO: Δiw_to_Cτ

include("../../../src/includes.jl")

using DelimitedFiles, JSON, Serialization, Interpolations

ρ₀(ϵ, t=5.) = sqrt(4*t^2-ϵ^2) / (2*π*t^2)

const t = 5.

function J(D::Real, ω::Real)
	return sqrt(1-(ω/D)^2) / π
end

spectrum_func(D=10.) = SpectrumFunction(ω -> J(D, ω), lb = -D, ub = D)


function main(; β=10., δτ=0.1, U=1., ϵ_d=U/2, chi=100)
	N = round(Int, β / δτ)
	D = 10.
	μ = 0.
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	lattice = GrassmannLattice(N=N, δτ=β/N, bands=2, contour=:imag)

	# initial guess for Δiw
	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, U=U, μ=-ϵ_d)
	corr = correlationfunction(bath, lattice)
	Δiw = toulouse_Δiw(bath)

	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
	mpsKs = [mpsK]
	for band in 1:lattice.bands
		mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
	end

	for i in 1:10
		println("the $i-th DMFT iteration...")
		corr = Δiw_to_Cτ(Δiw, β=β, N=N)
		mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)

		# compute observable
		cache = environments(lattice, mpsKs, mpsI1, mpsI2)

		gτ = cached_Gτ(lattice, mpsKs, mpsI1, mpsI2, cache=cache)
		# Gτ to Giw
		interp=linear_interpolation(0:δτ:β, gτ)
		gτ′ = interp.(0:0.01:β)
		Giw = Gτ_to_Giw(gτ′, β=β)

		###########general case##############
		# compute Σiw from Giw and G₀iw
		# compute lattice Giw
		# compute lattice G₀iw

		###########bethe lattice#############
		ws = ifrequencies(Giw, β=β)
		G0iw = [1 / (im*w - ϵ_d - t^2 * Gj) for (w, Gj) in zip(ws, Giw)]
		Δiw_new = [im*w - ϵ_d - 1/G0j for (w, G0j) in zip(ws, G0iw)]

		err = norm(Δiw_new - Δiw) / norm(Δiw)
		println("error is ", err)
		
		Δiw = Δiw_new
	end

	return Δiw

end

