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
			k = lattice.k-1
			mpsa =	vacuumstate(lattice)
			for fi in (:+, :-)
				η₁, η₂ = branch(corr, fi, :+), branch(corr, fi, :-)
				for i in 1:k, band in 1:lattice.bands
					mps1 = vacuumstate(lattice)
					for j in 1:k, band2 in 1:lattice.bands
						pos1a, pos1b = index(lattice, i+1, conj=true, band=band, branch=fi), index(lattice, i, conj=false, band=band, branch=fi)
						for fj in (:+, :-)
							η = (fj == :+) ? η₁ : η₂
							if (i == j) && (band == band2) && (fi == fj)
								t = exp(GTerm(pos1a, pos1b, coeff=η[i, j]))
							else
								pos2a, pos2b = index(lattice, j+1, conj=true, band=band2, branch=fj), index(lattice, j, conj=false, band=band2, branch=fj)
								t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=η[i, j]))
							end
							apply!(t, mps1)
							canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))	
						end
					end

					mult!(mpsa, mps1, trunc=trunc)
				end						
			end
			mpsb = retardedinteractdynamics(lattice, corr, trunc=trunc)
			@test distance(mpsa, mpsb) / norm(mpsa) <= rtol
		end
	end

	# 2 bands

	for N in (2, 3)
		for ordering in real_grassmann_orderings
			lattice = GrassmannLattice(N=N, δt=0.05, bands=2, contour=:real, order=1, ordering=ordering)
			corr = fermionic_Ct(f, β=β, N=lattice.N, t=lattice.t)
			k = lattice.k-1
			mpsa =	vacuumstate(lattice)

			# IF 11 and 22
			for band in 1:lattice.bands
				for fi in (:+, :-)
					η₁, η₂ = branch(corr, fi, :+), branch(corr, fi, :-)
					for i in 1:k, band in 1:lattice.bands
						mps1 = vacuumstate(lattice)
						for j in 1:k, band2 in 1:lattice.bands
							pos1a, pos1b = index(lattice, i+1, conj=true, band=band, branch=fi), index(lattice, i, conj=false, band=band, branch=fi)
							for fj in (:+, :-)
								η = (fj == :+) ? η₁ : η₂
								if (i == j) && (band == band2) && (fi == fj)
									t = exp(GTerm(pos1a, pos1b, coeff=η[i, j]))
								else
									pos2a, pos2b = index(lattice, j+1, conj=true, band=band2, branch=fj), index(lattice, j, conj=false, band=band2, branch=fj)
									t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=η[i, j]))
								end
								apply!(t, mps1)
								canonicalize!(mps1, alg=Orthogonalize(trunc = trunc))	
							end
						end

						mult!(mpsa, mps1, trunc=trunc)
					end						
				end
			end
			# IF 12
			for i in 1:k, fi in (:+, :-) 
				for j in 1:k, fj in (:+, :-)
					pos1a, pos1b = index(lattice, i+1, conj=true, band=1, branch=fi), index(lattice, i, conj=false, band=1, branch=fi)
					pos2a, pos2b = index(lattice, j+1, conj=true, band=2, branch=fj), index(lattice, j, conj=false, band=2, branch=fj)
					coef = index(corr, i, j, b1=fi, b2=fj) + index(corr, j, i, b1=fi, b2=fj)
					t = exp(GTerm(pos1a, pos1b, pos2a, pos2b, coeff=coef))
					apply!(t, mpsa)
					canonicalize!(mpsa, alg=Orthogonalize(TK.SVD(), trunc))
		
				end			
			end			

			mpsb = retardedinteractdynamics(lattice, corr, trunc=trunc)
			@test distance(mpsa, mpsb) / norm(mpsa) <= rtol
		end
	end	

end
