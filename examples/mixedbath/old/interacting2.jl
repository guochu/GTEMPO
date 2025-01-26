# push!(LOAD_PATH, "../../../src")
# using GTEMPO


include("../../../src/includes.jl")
using DelimitedFiles, JSON, Serialization


spectrum_func() = Leggett(d=3, ωc=1, α=0.4)

# spectrum_func() = DiracDelta(ω₀=1, α=0.5)

function main_analytic(U, ϵ_d=U/2)
	# ϵ_d = 0.5
	δτ=0.01
	N = 20
	β = N * δτ

	return independentbosons_Gτ(spectrum_func(), β=β, ϵ_d=-ϵ_d, U=U, N=N, nbands=2)
end

# function retardedinteractdynamics2!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
# 	corr = corr1.data
# 	k = lattice.k-1
# 	for i in 1:k, band1 in 1:lattice.bands
# 		pos1a, pos1b = index(lattice, i+1, conj=true, band=band1), index(lattice, i, conj=false, band=band1)
# 		for j in 1:k, band2 in 1:lattice.bands
# 			if (i == j) && (band1 == band2)
# 				coef = corr[i, j]
# 				t = exp(GTerm(pos1a, pos1b, coeff=coef))
# 			else
# 				pos2a, pos2b = index(lattice, j+1, conj=true, band=band2), index(lattice, j, conj=false, band=band2)
# 				coef = corr[i, j]
# 				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
# 			end
# 			apply!(t, gmps)
# 			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
# 		end
# 	end
# 	return gmps	
# end

# function retardedinteractdynamics2!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
# 	@assert lattice.bands == 2
# 	corr = corr1.data
# 	k = lattice.k-1
# 	for i in 1:k, band1 in 1:lattice.bands
# 		pos1a, pos1b = index(lattice, i+1, conj=true, band=band1), index(lattice, i, conj=false, band=band1)
# 		for j in 1:k, band2 in band1:lattice.bands
# 			if (i == j) && (band1 == band2)
# 				coef = corr[i, j]
# 				t = exp(GTerm(pos1a, pos1b, coeff=coef))
# 			else
# 				pos2a, pos2b = index(lattice, j+1, conj=true, band=band2), index(lattice, j, conj=false, band=band2)
# 				# coef = 2 * corr[i, j]
# 				coef = corr[i, j] + corr[j, i] 
# 				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
# 			end
# 			apply!(t, gmps)
# 			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
# 		end
# 	end
# 	return gmps	
# end

function retardedinteractdynamics2!(gmps::GrassmannMPS, lattice::ImagGrassmannLattice1Order, corr1::ImagCorrelationFunction; trunc::TruncationScheme=DefaultITruncation)
	@assert lattice.bands == 2
	corr = corr1.data
	k = lattice.k-1
	for i in 1:k, j in 1:k
		for band in 1:lattice.bands
			pos1a, pos1b = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
			pos2a, pos2b = index(lattice, j+1, conj=true, band=band), index(lattice, j, conj=false, band=band)
			if i == j
				coef = corr[i, j]
				t = exp(GTerm(pos1a, pos1b, coeff=coef))
			else
				coef = corr[i, j]
				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
			end
			apply!(t, gmps)
			canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
		end
	end

	for i in 1:k, j in 1:k
		pos1a, pos1b = index(lattice, i+1, conj=true, band=1), index(lattice, i, conj=false, band=1)
		pos2a, pos2b = index(lattice, j+1, conj=true, band=2), index(lattice, j, conj=false, band=2)
		coef = corr[i, j] + corr[j, i]
		t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
		apply!(t, gmps)
		canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))

		# pos1a, pos1b = index(lattice, i+1, conj=true, band=2), index(lattice, i, conj=false, band=2)
		# pos2a, pos2b = index(lattice, j+1, conj=true, band=1), index(lattice, j, conj=false, band=1)
		# coef = corr[i, j] 
		# t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
		# apply!(t, gmps)
		# canonicalize!(gmps, alg=Orthogonalize(TK.SVD(), trunc))
	end
	return gmps	
end

retardedinteractdynamics2(lattice::ImagGrassmannLattice1Order, corr::ImagCorrelationFunction; kwargs...) = retardedinteractdynamics2!(vacuumstate(lattice), lattice, corr; kwargs...)

function main(U, ϵ_d=U/2)
	# ϵ_d = 0.5
	δτ=0.01
	N = 20
	β = N * δτ
	chi = 120

	trunc = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)
	truncK = truncdimcutoff(D=chi, ϵ=1.0e-12, add_back=0)


	
	lattice = GrassmannLattice(N=N, δτ=δτ, contour=:imag, order=1, bands=2)

	bath = bosonicbath(spectrum_func(), β=β)

	corr = correlationfunction(bath, lattice)
	# println(corr.data.ηⱼₖ, " ", corr.data.ηₖⱼ)
	g₀ = exp(-δτ*ϵ_d)
	# η1, η2 = corr.data.ηⱼₖ, corr.data.ηₖⱼ
	# η1[1] = η1[1] + η2[1]
	# η2[1] = 0
	# η1 = g₀^2 .* (exp.(η1) .- 1)
	# η2 = g₀^2 .* (exp.(η2) .- 1)
	# corr = ImagCorrelationFunction(CorrelationMatrix(η1, η2))


	println("computing MPS-IF...")
	# @time mpsI1 = retardedinteractdynamics(lattice, corr, trunc=trunc, band=1)
	@time mpsI = retardedinteractdynamics2(lattice, corr, trunc=trunc)

	println("bond dimension of mpsI is ", bond_dimension(mpsI))

	fbath = fermionicbath(semicircular(), β=β, μ=0)
	exact_model = AndersonIM(U=U, μ=-ϵ_d)
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