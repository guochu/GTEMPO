# push!(LOAD_PATH, "../../../src")
# using GTEMPO


include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


spectrum_func() = Leggett(d=3, ωc=10)

function main_analytic(ϵ_d)
	# ϵ_d = 0.5
	δτ=0.05
	N = 10
	β = N * δτ

	return independentbosons_Gτ(spectrum_func(), β=β, ϵ_d=-ϵ_d, N=N)

end

function main(ϵ_d)
	# ϵ_d = 0.5
	δτ=0.05
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
	exact_model = SISB(fbath, U=0., μ=-ϵ_d)
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