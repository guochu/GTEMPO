# push!(LOAD_PATH, "../../../src")
# using GTEMPO


include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


# spectrum_func() = Leggett(d=3, ωc=1)

spectrum_func() = DiracDelta(ω₀=1, α=0.5)

function main_imag_analytic(ϵ_d)
	# ϵ_d = 0.5
	δτ=0.01
	N = 10
	β = N * δτ

	return independentbosons_Gτ(spectrum_func(), β=β, ϵ_d=-ϵ_d, N=N)
end

function main_imag(ϵ_d)
	# ϵ_d = 0.5
	δτ=0.01
	N = 10
	β = N * δτ
	chi = 100

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)


	
	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1)

	bath = bosonicbath(spectrum_func(), β=β)
	corr = correlationfunction(bath, lattice)

	println("computing MPS-IF...")
	@time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	println("bond dimension of mpsI is ", bond_dimension(mpsI))

	fbath = fermionicbath(semicircular(), β=β, μ=0)
	exact_model = AndersonIM(fbath, U=0., μ=-ϵ_d)
	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))
	mpsKs = [mpsK]
	for band in 1:lattice.bands
		mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
	end
	println("bond dimension of mpsK is ", bond_dimension(mpsKs[1]), ", number of Ks ", length(mpsKs))




	cache = environments(lattice, mpsKs, mpsI)
	@time g = cached_Gτ(lattice, mpsKs, mpsI, cache=cache)

	return g
end

function main_real_analytic(ϵ_d)
	# ϵ_d = 0.5
	β = 0.1
	δt=0.01
	N = 10
	t = N * δt
	g1 = [independentbosons_greater(spectrum_func(), tj, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
	g2 = [independentbosons_lesser(spectrum_func(), tj, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
	return g1, g2
end

function main_real(ϵ_d)
	# ϵ_d = 0.5
	β = 0.1
	δt=0.01
	Nt = 10
	t = Nt * δt
	chi = 100

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(N=Nt, δt=δt, contour=:real, order=1)

	bath = bosonicbath(spectrum_func(), β=β)
	corr = correlationfunction(bath, lattice)

	println("computing MPS-IF...")
	@time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	println("bond dimension of mpsI is ", bond_dimension(mpsI))

	fbath = fermionicbath(semicircular(), β=β, μ=0)
	exact_model = AndersonIM(fbath, U=0., μ=-ϵ_d)
	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	mpsK = systhermalstate!(mpsK, lattice, exact_model, trunc=trunc, δτ=0.001)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))




	cache = environments(lattice, mpsK, mpsI)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time g₁ = [cached_greater(lattice, k, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	@time g₂ = [cached_lesser(lattice, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]

	return -im*g₁, im*g₂
end

function main_mixed(ϵ_d)
	# ϵ_d = 0.5
	β = 0.1
	δτ = 0.01
	Nτ = round(Int, β/δτ)
	δt=0.01
	Nt = 10
	t = Nt * δt
	chi = 120

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)

	
	lattice = GrassmannLattice(Nt=Nt, δt=δt, Nτ=Nτ, δτ=δτ, contour=:mixed, order=1)

	bath = bosonicbath(spectrum_func(), β=β)
	corr = correlationfunction(bath, lattice)

	println("computing MPS-IF...")
	@time mpsI = retardedinteractdynamics(lattice, corr, trunc=trunc)

	println("bond dimension of mpsI is ", bond_dimension(mpsI))

	fbath = fermionicbath(semicircular(), β=β, μ=0)
	exact_model = AndersonIM(fbath, U=0., μ=-ϵ_d)
	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))
	for band in 1:lattice.bands
		mpsK = boundarycondition!(mpsK, lattice, band=band)
	end
	# mpsK = systhermalstate!(mpsK, lattice, exact_model, trunc=trunc, δτ=0.001)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))




	cache = environments(lattice, mpsK, mpsI)
	# @time gt = [cached_greater(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# @time lt = [cached_lesser(lattice, j, mpsK, mpsI, cache=cache) for j in 1:lattice.kt]
	# return gt, lt

	@time g₁ = [cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:+, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	@time g₂ = [cached_Gm(lattice, 1, k, mpsK, mpsI, c1=true, c2=false, b1=:-, b2=:+, band=1, cache=cache) for k in 1:Nt+1]
	@time g₃ = [cached_Gm(lattice, k, 1, mpsK, mpsI, c1=false, c2=true, b1=:τ, b2=:τ, band=1, cache=cache) for k in 1:Nτ+1]

	return -im*g₁, im*g₂, real(g₃)
end



