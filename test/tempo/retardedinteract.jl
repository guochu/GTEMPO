println("------------------------------------")
println("|    PartialIF-Retarded Interact   |")
println("------------------------------------")


@testset "Retarded Interact IF: imaginary time" begin
	trunc = truncdimcutoff(D=300, ϵ=1.0e-8, add_back=0)
	rtol = 1.0e-5
	# 1 band
	for N in (2,3)
		for ordering in imag_grassmann_orderings

			# println("bands = ", bands, ", N = ", N)
			lattice = GrassmannLattice(N=N, δτ=0.05, bands=1, contour=:imag, order=1, ordering=ordering)
			η = randn(Float64, lattice.k-1, lattice.k-1)
			band = 1
			mpsa = vacuumstate(lattice)
			for i in 1:lattice.k-1
				mps1 = vacuumstate(lattice)
				for j in 1:lattice.k-1
					if (i == j) 
						pos1, pos2 = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
						t = exp(GTerm(pos1, pos2, coeff=η[i, j]))
					else
						pos1a, pos1b = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
						pos2a, pos2b = index(lattice, j+1, conj=true, band=band), index(lattice, j, conj=false, band=band)
						t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=η[i, j]))
					end
					mps1 = t * mps1
					canonicalize!(mps1, alg=Orthogonalize(TK.SVD(), trunc))
				end	
				mult!(mpsa, mps1, trunc=trunc)
			end
			mpsb = retardedinteractdynamics(lattice, ImagCorrelationFunction(η), trunc=trunc)
			@test distance(mpsa, mpsb) / norm(mpsa) <= rtol

			mpsc = retardedinteractdynamics_naive(lattice, ImagCorrelationFunction(η), trunc=trunc)
			@test distance(mpsa, mpsc) / norm(mpsa) <= rtol
		end
	end		


	# 2 band
	for N in (2,3)
		for ordering in imag_grassmann_orderings
			# println("ordering is ", ordering, ", N = ", N)
			# println("bands = ", bands, ", N = ", N)
			lattice = GrassmannLattice(N=N, δτ=0.05, bands=2, contour=:imag, order=1, ordering=ordering)
			η = randn(Float64, lattice.k-1, lattice.k-1)

			mpsa =	vacuumstate(lattice)
			k = lattice.k-1
			for i in 1:k, j in 1:k
				for band in 1:lattice.bands
					pos1a, pos1b = index(lattice, i+1, conj=true, band=band), index(lattice, i, conj=false, band=band)
					pos2a, pos2b = index(lattice, j+1, conj=true, band=band), index(lattice, j, conj=false, band=band)
					if i == j
						coef = η[i, j]
						t = exp(GTerm(pos1a, pos1b, coeff=coef))
					else
						coef = η[i, j]
						t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
					end
					apply!(t, mpsa)
					canonicalize!(mpsa, alg=Orthogonalize(TK.SVD(), trunc))
				end
			end

			for i in 1:k, j in 1:k
				pos1a, pos1b = index(lattice, i+1, conj=true, band=1), index(lattice, i, conj=false, band=1)
				pos2a, pos2b = index(lattice, j+1, conj=true, band=2), index(lattice, j, conj=false, band=2)
				coef = η[i, j] + η[j, i]
				t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
				apply!(t, mpsa)
				canonicalize!(mpsa, alg=Orthogonalize(TK.SVD(), trunc))
			end

			mpsb = retardedinteractdynamics(lattice, ImagCorrelationFunction(η), trunc=trunc)
			@test distance(mpsa, mpsb) / norm(mpsa) <= rtol

			mpsc = retardedinteractdynamics_naive(lattice, ImagCorrelationFunction(η), trunc=trunc)
			@test distance(mpsa, mpsc) / norm(mpsa) <= rtol
		end
	end	

end


@testset "Retarded Interact IF: real time" begin
	D = 1
	β=0.2
	J(ε) = D/(ε^2+D^2)/pi
	f = SpectrumFunction(J)
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)
	rtol = 1.0e-5

	# 1 band
	for N in (2, 3, 4)
		for ordering in real_grassmann_orderings
			# println("bands=", bands, ", N=", N)
			lattice = GrassmannLattice(N=N, δt=0.05, bands=1, contour=:real, order=1, ordering=ordering)
			corr = fermionic_Ct(f, β=β, N=lattice.N, t=lattice.t)
			mpsa = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
			mpsb = retardedinteractdynamics(lattice, corr, trunc=trunc)
			@test distance(mpsa, mpsb) / norm(mpsa) <= rtol
		end
	end

	# 2 bands

	for N in (2, 3)
		for ordering in real_grassmann_orderings
			lattice = GrassmannLattice(N=N, δt=0.05, bands=2, contour=:real, order=1, ordering=ordering)
			corr = fermionic_Ct(f, β=β, N=lattice.N, t=lattice.t)

			mpsa = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)	
			mpsb = retardedinteractdynamics(lattice, corr, trunc=trunc)
			@test distance(mpsa, mpsb) / norm(mpsa) <= rtol
		end
	end	

end

@testset "Retarded Interact IF: mixed time" begin
	D = 1
	β=0.2
	Nτ = 2
	δτ = β / Nτ
	J(ω) = (D/(2*pi)) * sqrt(1 - (ω/D)^2 ) * 0.1
	f = SpectrumFunction(J, lb = -D, ub = D)
	trunc = truncdimcutoff(D=300, ϵ=1.0e-6, add_back=0)
	rtol = 1.0e-5

	for bands in (1, 2)
		for N in (1, 2, 3)
			for ordering in mixed_grassmann_orderings
				lattice = GrassmannLattice(Nt=N, δt=0.05, δτ=δτ, Nτ=Nτ, bands=bands, contour=:mixed, order=1, ordering=ordering)
				corr = fermionic_Cm(f, β=β, Nt=lattice.Nt, t=lattice.t, Nτ=lattice.Nτ)

				mpsa = retardedinteractdynamics_naive(lattice, corr, trunc=trunc)
				mpsb = retardedinteractdynamics(lattice, corr, trunc=trunc)
				@test distance(mpsa, mpsb) / norm(mpsa) <= rtol
			end
		end
	end

end


