# push!(LOAD_PATH, "../../../src")
# using GTEMPO


include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


spectrum_func(α) = Leggett(d=3, ωc=10, α=α)

function main_analytic(U, ϵ_d=U/2)
	# ϵ_d = 0.5
	δτ=0.05
	N = 20
	β = N * δτ

	return independentbosons_Gτ(spectrum_func(1), β=β, ϵ_d=-ϵ_d, U=U, N=N, nbands=2)
end

function partialif_retardedinteract2(lattice::ImagGrassmannLattice1Order, i::Int, cols::AbstractVector, g; band::Int=1, trunc::TruncationScheme=DefaultMPOTruncation)
	@assert lattice.k - 1 == length(cols)

	mps1 = vacuumstate(lattice)
	for j in 1:lattice.k-1, band2 in 1:lattice.bands
		if (i == j) && (band == band2)
			pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
			coef = cols[j]
			# coef = g^2 * (exp(cols[j])-1)
			t = exp(GTerm(pos1, pos2, coeff=coef))
		else
			pos1a, pos1b = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
			pos2a, pos2b = index(lattice, j+1, conj=true, band=band2), index(lattice, j, conj=false, band=band2)
			coef = g^2 * (exp(cols[j])-1)
			t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
		end
		apply!(t, mps1)
		canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))
	end	

	return mps1
end

function retardedinteractdynamics2!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction, g; band::Int=1, trunc::TruncationScheme=DefaultITruncation)
	corr = corr1.data
	k = lattice.k-1
	for i in 1:k, band2 in 1:lattice.bands
		tmp = partialif_retardedinteract2(lattice, i, view(corr, i, 1:k), g, band=band)
		gmps = mult!(gmps, tmp, trunc=trunc)
	end
	return gmps	
end
retardedinteractdynamics2(lattice::ImagGrassmannLattice1Order, corr::ImagCorrelationFunction, g; kwargs...) = retardedinteractdynamics2!(vacuumstate(lattice), lattice, corr, g; kwargs...)

function main(U, ϵ_d=U/2)
	# ϵ_d = 0.5
	δτ=0.05
	N = 20
	β = N * δτ
	chi = 120

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-10, add_back=0)


	
	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=2)

	bath = bosonicbath(spectrum_func(1), β=β)

	corr = correlationfunction(bath, lattice)
	g₀ = exp(-δτ*ϵ_d)
	# η1, η2 = corr.data.ηⱼₖ, corr.data.ηₖⱼ
	# η1[1] = η1[1] + η2[1]
	# η2[1] = 0
	# η1 = g₀^2 .* (exp.(η1) .- 1)
	# η2 = g₀^2 .* (exp.(η2) .- 1)
	# corr = ImagCorrelationFunction(CorrelationMatrix(η1, η2))


	println("computing MPS-IF...")
	# @time mpsI1 = retardedinteractdynamics(lattice, corr, trunc=trunc, band=1)
	@time mpsI1 = retardedinteractdynamics2(lattice, corr, g₀, trunc=trunc, band=1)
	@time mpsI2 = swapband(mpsI1, lattice, 1, 2, trunc=trunc)

	println("bond dimension of mpsI is ", bond_dimension(mpsI1), " ", bond_dimension(mpsI2))

	fbath = fermionicbath(semicircular(), β=β, μ=0)
	exact_model = SISB(fbath, U=U, μ=-ϵ_d)
	@time mpsK = sysdynamics(lattice, exact_model, trunc=truncK)
	println("bond dimension of mpsK is ", bond_dimension(mpsK))
	println("mpsK scale is ", scaling(mpsK))
	mpsKs = [mpsK]
	for band in 1:lattice.bands
		mpsKs = boundarycondition_branching(mpsKs, lattice, band=band)
	end
	println("bond dimension of mpsK is ", bond_dimension(mpsKs[1]), ", number of Ks ", length(mpsKs))




	cache = environments(lattice, mpsKs, mpsI1, mpsI2)
	@time g = cached_Gτ(lattice, mpsKs, mpsI1, mpsI2, cache=cache)

	return g
end