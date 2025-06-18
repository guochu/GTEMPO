println("------------------------------------")
println("|    PartialIF-Retarded Interact   |")
println("------------------------------------")


@testset "Retarded Interact IF: imaginary time" begin
	trunc = truncdimcutoff(D=300, ϵ=1.0e-8, add_back=0)
	rtol = 1.0e-5
	# 1 band
	for N in (2,3)
		# println("bands = ", bands, ", N = ", N)
		lattice = FockLattice(N=N, δτ=0.05, bands=1, contour=:imag, order=1)
		η = randn(Float64, lattice.N, lattice.N)
		band = 1
		mpsa = vacuumstate(lattice)
		for i in 1:lattice.N
			mps1 = vacuumstate(lattice)
			for j in 1:lattice.N
				coef = exp(η[i, j]) 
				if (i == j) 
					pos1 = index(lattice, i, band=band)
					t = ExpNTerm(pos1, coeff=coef)
				else
					pos1, pos2 = index(lattice, i, band=band), index(lattice, j, band=band)
					t = ExpNTerm(pos1, pos2, coeff=coef)
				end
				apply!(t, mps1)
				canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))
			end	
			mult!(mpsa, mps1, trunc=trunc)
		end

		mpsb = hybriddynamics_naive(lattice, ImagCorrelationFunction(η), trunc=trunc)
		@test distance(mpsa, mpsb) / norm(mpsa) <= rtol

		mpsc = hybriddynamics(lattice, ImagCorrelationFunction(η), trunc=trunc)
		@test distance(mpsa, mpsc) / norm(mpsa) <= rtol
	end		


	# 2 band
	for N in (2,3)
		# println("ordering is ", ordering, ", N = ", N)
		# println("bands = ", bands, ", N = ", N)
		lattice = FockLattice(N=N, δτ=0.05, bands=2, contour=:imag, order=1)
		η = randn(Float64, lattice.N, lattice.N)

		mpsa =	vacuumstate(lattice)
		k = lattice.N
		for i in 1:k, j in 1:k
			for band in 1:lattice.bands
				pos1, pos2 = index(lattice, i, band=band), index(lattice, j, band=band)
				coef = exp(η[i, j]) 
				if pos1 == pos2
					t = ExpNTerm(pos1, coeff=coef)
				else
					t = ExpNTerm(pos1, pos2, coeff=coef)
				end
				apply!(t, mpsa)
				canonicalize!(mpsa, alg=Orthogonalize(TK.SVD(), trunc))
			end
		end

		for i in 1:k, j in 1:k
			pos1, pos2 = index(lattice, i, band=1), index(lattice, j, band=2)
			coef = exp(η[i, j] + η[j, i]) 
			t = ExpNTerm(pos1, pos2, coeff=coef)
			apply!(t, mpsa)
			canonicalize!(mpsa, alg=Orthogonalize(TK.SVD(), trunc))
		end

		mpsb = hybriddynamics_naive(lattice, ImagCorrelationFunction(η), trunc=trunc)
		@test distance(mpsa, mpsb) / norm(mpsa) <= rtol

		mpsc = hybriddynamics(lattice, ImagCorrelationFunction(η), trunc=trunc)
		@test distance(mpsa, mpsc) / norm(mpsa) <= rtol
	end	

end

@testset "Retarded Interact IF: real time" begin
	D = 1
	β=0.2
	f = Leggett(d=3, ωc=5, α=0.4)
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)
	rtol = 1.0e-5

	# 1 band
	for N in (2, 3, 4)
		# println("bands=", bands, ", N=", N)
		lattice = FockLattice(N=N, δt=0.05, bands=1, contour=:real, order=1)
		corr = bosonic_Δt(f, β=β, N=lattice.N, t=lattice.t)
		mpsa = hybriddynamics_naive(lattice, corr, trunc=trunc)
		mpsb = hybriddynamics(lattice, corr, trunc=trunc)
		@test distance(mpsa, mpsb) / norm(mpsa) <= rtol
	end

	# 2 bands

	for N in (2, 3)
		lattice = FockLattice(N=N, δt=0.05, bands=2, contour=:real, order=1)
		corr = bosonic_Δt(f, β=β, N=lattice.N, t=lattice.t)

		mpsa = hybriddynamics_naive(lattice, corr, trunc=trunc)	
		mpsb = hybriddynamics(lattice, corr, trunc=trunc)
		@test distance(mpsa, mpsb) / norm(mpsa) <= rtol
	end	
end

@testset "Retarded Interact IF: mixed time" begin
	D = 1
	β=0.2
	Nτ = 2
	δτ = β / Nτ
	f = Leggett(d=1, ωc=5, α=0.4)
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)
	rtol = 1.0e-4

	for bands in (1, 2)
		for N in (1, 2, 3)
			# println("N=", N, ", bands=", bands)
			lattice = FockLattice(Nt=N, δt=0.05, δτ=δτ, Nτ=Nτ, bands=bands, contour=:mixed, order=1)
			corr = bosonic_Δm(f, β=β, Nt=lattice.Nt, t=lattice.t, Nτ=lattice.Nτ)

			mpsa = hybriddynamics_naive(lattice, corr, trunc=trunc)
			mpsb = hybriddynamics(lattice, corr, trunc=trunc)
			@test distance(mpsa, mpsb) / norm(mpsa) <= rtol
		end
	end

end
