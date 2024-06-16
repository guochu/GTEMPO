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

	mpath = "beta$(round(Int, β))t$(round(Int, t))/"

	lattice = GrassmannLattice(Nτ=Nτ, δτ=δτ, δt=δt, Nt=Nt, bands=2, contour=:mixed)

	# initial guess for Δiw
	bath = fermionicbath(spectrum_func(D), β=β, μ=0)
	exact_model = SISB(bath, U=U, μ=-ϵ_d)

	# initial guess for Δiw
	lb = -5.
	ub = 5.
	dw = 1.0e-4
	freqs = collect(frequencies(lb=lb, ub=ub, dw=dw))
	Jw =  [toulouse_Jw(spectrum_func(D), ω) for ω in freqs]


	mpsK = sysdynamics(lattice, exact_model, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition(mpsK, lattice, band=band)
	end

	for i in 1:10
		println("the $i-th DMFT iteration...")

		spec = SpectrumFunction(freqs, Jw)
		bath = fermionicbath(spec, β=β, μ=0)
		corr = correlationfunction(bath, lattice)


		mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)

		# compute observable
		cache = environments(lattice, mpsK, mpsI1, mpsI2)

		gt = [cached_greater(lattice, k, mpsK, mpsI1, mpsI2, cache=cache) for k in 1:lattice.kt]
		lt = [cached_lesser(lattice, k, mpsK, mpsI1, mpsI2, cache=cache) for k in 1:lattice.kt]

		retarded = -im * gt - im * lt
		# println(retarded)
		# linear prediction
		retarded′ = linear_predict(retarded, δt, tol=1.0e-9)
		println(retarded′[end-5:end])
		t_all = δt * (length(retarded′)-1)
		println("tall is ", t_all)

		# Gt to Gw
		interp=linear_interpolation(0:δt:t_all, retarded′)
		δt′ = 0.001
		gf′ = interp.(0:δt′:t_all)
		Gw = Gt_to_Gw(gf′, δt′, lb=lb, ub=ub, dw=dw)
		println("Gw[1]=", Gw[1], ", Gw[end]=", Gw[end])

		###########general case##############
		# compute Σiw from Giw and G₀iw
		# compute lattice Giw
		# compute lattice G₀iw

		###########bethe lattice#############
		Δw_new = [Dh^2 * Gj for (w, Gj) in zip(freqs, Gw)]
		Jw_new = Gw_to_Aw(Δw_new, verbosity=0)

		err = norm(Jw_new - Jw) / norm(Jw)
		println("error is ", err)

		writedlm(mpath*"Gt-$i.dat", hcat(0:δt′:t_all, real(gf′), imag(gf′)))
		writedlm(mpath*"Giw-$i.dat", hcat(freqs, real(Gw), imag(Gw)))
		
		Jw = Jw_new

		Aw = Gw_to_Aw(Gw, verbosity=0)
		nrm = sum(Aw) * dw
		println("sum(Aw)=", nrm)
		Aw ./= nrm

		Gtau = Aw_to_Gτ(Aw, β=β, lb=lb, ub=ub, dw=dw, δτ=δτ)
		writedlm(mpath*"Aw-$i.dat", hcat(freqs, Aw))
		writedlm(mpath*"Gtau-$i.dat", hcat(0:δτ:β, Gtau))
	end

	return Jw

end

