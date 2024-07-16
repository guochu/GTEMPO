# push!(LOAD_PATH, "../../../src")

# using GTEMPO
# using GTEMPO: Δiw_to_Cτ

include("../../../src/includes.jl")

using DelimitedFiles, JSON, Serialization, Interpolations

ρ₀(ϵ, D=1.) = sqrt(1-(ϵ/D)^2) * (D/π)


spectrum_func(D=1.) = SpectrumFunction(ω -> 0.1 * ρ₀(ω, D), lb = -D, ub = D)


function quench_sysdynamics(lattice, bath; ϵ_d, U, trunc)
	model1 = SISB(bath, U=U, μ=-ϵ_d)
	model2 = SISB(bath, U=0, μ=-ϵ_d)
	mpsK = sysdynamics_forward!(vacuumstate(lattice), lattice, model1, trunc=trunc)
	mpsK = sysdynamics_backward!(mpsK, lattice, model1, trunc=trunc)
	mpsK = sysdynamics_imaginary!(mpsK, lattice, model2, trunc=trunc)
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	return mpsK
end


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


	mpsK = quench_sysdynamics(lattice, bath, ϵ_d=ϵ_d, U=U, trunc=trunc)

	spec = SpectrumFunction(freqs, Jw)
	bath = fermionicbath(spec, β=β, μ=0)
	corr = correlationfunction(bath, lattice)

	for i in 1:10
		println("the $i-th DMFT iteration...")

		mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
		mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)

		# compute observable
		cache = environments(lattice, mpsK, mpsI1, mpsI2)

		@time gfs = cached_Gm(lattice, mpsK, mpsI1, mpsI2, cache=cache) 

		###########general case##############
		# compute Σiw from Giw and G₀iw
		# compute lattice Giw
		# compute lattice G₀iw

		###########bethe lattice#############
		Δw_new = [Dh^2 * Gj for Gj in gfs]
		for i in 1:2, j in 1:2
			Δw_new[i, j] .*= (-im*δt^2)
		end
		Δw_new[3,3] .*= (-im*δτ^2)
		for i in 1:2
			Δw_new[i, 3] .*= (-im*δt * δτ)
			Δw_new[3, i] .*= (-im*δt * δτ)
		end

		corr_new = Cm(Δw_new)

		err = [norm(branch(corr_new, b1=b1, b2=b2) - branch(corr, b1=b1, b2=b2)) / norm(branch(corr, b1=b1, b2=b2)) for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)]
		println("error is ", err)
		
		corr = corr_new
	end

	return corr

end

function test(; β=10., δτ=0.1, t=2., δt=0.05, U=1., ϵ_d=U/2, chi=100)
	Nτ = round(Int, β / δτ)
	Nt = round(Int, t / δt)
	total_t = 2.
	δt = 0.05
	
	Dh = 0.5
	D = 2 * Dh
	μ = 0.

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	mpath = "beta$(round(Int, β))t$(round(Int, t))/"

	lattice = GrassmannLattice(Nτ=Nτ, δτ=δτ, δt=δt, Nt=Nt, bands=1, contour=:mixed)

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



	spec = SpectrumFunction(freqs, Jw)
	bath = fermionicbath(spec, β=β, μ=0)
	corr = correlationfunction(bath, lattice)


	mpsI = hybriddynamics(lattice, corr, trunc=trunc, band=1)

	# compute observable
	cache = environments(lattice, mpsK, mpsI)

	@time gfs = cached_Gm(lattice, mpsK, mpsI, cache=cache) 

	return gfs
end


function test2(; β=10., δτ=0.1, t=2., δt=0.05, U=1., ϵ_d=U/2, chi=100)
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


	spec = SpectrumFunction(freqs, Jw)
	bath = fermionicbath(spec, β=β, μ=0)
	corr = correlationfunction(bath, lattice)


	mpsI1 = hybriddynamics(lattice, corr, trunc=trunc, band=1)
	mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)

	# compute observable
	cache = environments(lattice, mpsK, mpsI1, mpsI2)

	@time gfs = cached_Gm(lattice, mpsK, mpsI1, mpsI2, cache=cache) 

	###########general case##############
	# compute Σiw from Giw and G₀iw
	# compute lattice Giw
	# compute lattice G₀iw

	###########bethe lattice#############
	Δw_new = [Dh^2 * Gj for Gj in gfs]

	bhs = (:+, :-, :τ)
	for (i, b1) in enumerate(bhs)
		for (j, b2) in enumerate(bhs)
			corrj = branch(corr, b1=b1, b2=b2)
			println("error for btanch $i, $j is ", norm(Δw_new[i, j] - corrj) / norm(corrj) )
		end
	end
end