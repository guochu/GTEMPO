# push!(LOAD_PATH, "../../../src")

# using GTEMPO
# using GTEMPO: Δiw_to_Cτ

include("../../../src/includes.jl")

using DelimitedFiles, JSON, Serialization, Interpolations

ρ₀(ϵ, D=1.) = sqrt(1-(ϵ/D)^2) * (D/π)


spectrum_func(D=1.) = SpectrumFunction(ω -> ρ₀(ω, D), lb = -D, ub = D)


function main(; β=10., δτ=0.1, U=1., ϵ_d=U/2, chi=100)
	N = round(Int, β / δτ)
	t = 0.5
	D = 2 * t
	μ = 0.
	Nw = 1024
	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	mpath = "beta$(round(Int, β))/"

	lattice = GrassmannLattice(N=N, δτ=β/N, bands=2, contour=:imag)

	# initial guess for Δiw
	bath = fermionicbath(spectrum_func(2*t), β=β, μ=0)
	exact_model = AndersonIM(U=U, μ=-ϵ_d)

	# initial guess for Δiw
	Δiw =  toulouse_Δiw(bath, nmax=Nw) .* t^2

	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
	mpsKs = [mpsK]
	for band in 1:lattice.bands
		mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
	end

	for i in 1:10
		println("the $i-th DMFT iteration...")
		corr = Δiw_to_Cτ(Δiw, β=β, N=N)
		# bath2 = fermionicbath(SpectrumFunction(ω -> ρ₀(ω, D) * t^2, lb = -D, ub = D), β=β, μ=0)
		# corr2 = correlationfunction(bath2, lattice)
		# return corr, corr2

		mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)

		# compute observable
		cache = environments(lattice, mpsKs, mpsI1, mpsI2)

		gτ = -cached_Gτ(lattice, mpsKs, mpsI1, mpsI2, cache=cache)
		# Gτ to Giw
		interp=linear_interpolation(0:δτ:β, gτ)
		gτ′ = interp.(0:0.0001:β)
		Giw = Gτ_to_Giw(gτ′, β=β, nmax=Nw)

		###########general case##############
		# compute Σiw from Giw and G₀iw
		# compute lattice Giw
		# compute lattice G₀iw

		###########bethe lattice#############
		ws = ifrequencies(Giw, β=β)
		G0iw = [1 / (im*w + ϵ_d - t^2 * Gj) for (w, Gj) in zip(ws, Giw)]
		Δiw_new = [im*w + ϵ_d - 1/G0j for (w, G0j) in zip(ws, G0iw)]

		err = norm(Δiw_new - Δiw) / norm(Δiw)
		println("error is ", err)

		writedlm(mpath*"Gtau-$i.dat", gτ′)
		writedlm(mpath*"Giw-$i.dat", hcat(real(Giw), imag(Giw)))
		
		Δiw = Δiw_new
	end

	return Δiw

end

