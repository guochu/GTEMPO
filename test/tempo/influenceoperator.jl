println("------------------------------------")
println("|        InfluenceOperator         |")
println("------------------------------------")


@testset "InfluenceOperator: imaginary-time" begin
	β = 1.
	δτ = 0.2
	N = round(Int, β/δτ)
	bath = fermionicbath(spectrum_func(1), β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	algexpan = PronyExpansion(tol=1.0e-6)
	for bands in (1,2,3)
		for ordering in imag_grassmann_orderings
			if LayoutStyle(ordering) isa TimeLocalLayout

				lattice = GrassmannLattice(N=N, δτ=δτ, bands=bands, contour=:imag, ordering=ordering)
				corr = correlationfunction(bath, lattice)
				for band in 1:lattice.bands
					mpo1 = influenceoperator(lattice, corr, band=band, algexpan=algexpan)
					mps1 = mpo1 * vacuumstate(lattice)
					canonicalize!(mps1, alg=Orthogonalize(trunc=trunc))

					vc = vacuumstate(lattice)
					local mps2
					for i in 2:lattice.k
						for j in 2:lattice.k
							pos1, pos2 = index(lattice, i, conj=true, band=band), index(lattice, j, conj=false, band=band)
							if corr.data[i-1, j-1] != 0
								t = GTerm(pos1, pos2, coeff=corr.data[i-1, j-1])
								if @isdefined mps2
									mps2 += t * vc
								else
									mps2 = t * vc
								end
								canonicalize!(mps2, alg=Orthogonalize(trunc=trunc))
							end
						end
					end
					@test distance(mps1, mps2) / norm(mps2) < 1.0e-6

					dt = 0.01

					mps1 = dt * mps1 + vacuumstate(lattice)
					mpo2 = influenceoperatorexponential(lattice, corr, dt, WII(), band=band, algexpan=algexpan)
					mps0 = mpo2 * vacuumstate(lattice)

					@test distance(mps1, mps0) / norm(mps0) < dt

					for algmult in (SVDCompression(D=50), DMRG1(D=50), DMRG2(trunc=truncdimcutoff(D=50,ϵ=1.0e-6)))
						mps1 = differentialinfluencefunctional(lattice, corr, dt, WII(), algmult, band=band, algexpan=algexpan)
						_n = norm(mps1)
						mps2 = differentialinfluencefunctional(lattice, corr, dt, WI(), algmult, band=band, algexpan=algexpan)
						mps3 = differentialinfluencefunctional(lattice, corr, dt, ComplexStepper(WI()), algmult, band=band, algexpan=algexpan)
						mps4 = differentialinfluencefunctional(lattice, corr, dt, ComplexStepper(WII()), algmult, band=band, algexpan=algexpan)
						@test distance(mps1, mps0) / _n < dt
						@test distance(mps1, mps0) / _n < dt
						@test distance(mps1, mps0) / _n < dt
						@test distance(mps1, mps0) / _n < dt
					end

				end

			end

		end
	end
end


@testset "InfluenceOperator: real-time" begin
	β = 1.
	δt = 0.2
	N = round(Int, β/δt)
	bath = fermionicbath(spectrum_func2(1), β=β, μ=0)
	trunc = truncdimcutoff(D=50, ϵ=1.0e-6, add_back=0)
	algexpan = PronyExpansion(tol=1.0e-6)
	tol = 1.0e-5
	# only supports the following ordering currently
	ordering = 
	for ordering in [A1A1a1a1B1B1b1b1(), A1A1B1B1a1a1b1b1()]
		for bands in 1:3

			lattice = GrassmannLattice(N=N, δt=δt, bands=bands, contour=:real, ordering=ordering)
			corr = correlationfunction(bath, lattice)
			for band in 1:lattice.bands

				h_pp, h_pm, h_mp, h_mm = influenceoperator(lattice, corr, band=band, algexpan=algexpan)
				mps_pp = h_pp * vacuumstate(lattice)
				canonicalize!(mps_pp, alg=Orthogonalize(trunc=trunc))

				mps_pm = h_pm * vacuumstate(lattice)
				canonicalize!(mps_pm, alg=Orthogonalize(trunc=trunc))

				mps_mp = h_mp * vacuumstate(lattice)
				canonicalize!(mps_mp, alg=Orthogonalize(trunc=trunc))

				mps_mm = h_mm * vacuumstate(lattice)
				canonicalize!(mps_mm, alg=Orthogonalize(trunc=trunc))


				vc = vacuumstate(lattice)
				mps2_pp = nothing
				for i in 1:lattice.k
					for j in 1:lattice.k
						pos1, pos2 = index(lattice, i, branch=:+, conj=true, band=band), index(lattice, j, branch=:+, conj=false, band=band)
						if corr.G₊₊[i, j] != 0
							t = GTerm(pos1, pos2, coeff=corr.G₊₊[i, j])
							if !isnothing(mps2_pp)
								mps2_pp += t * vc
							else
								mps2_pp = t * vc
							end
							canonicalize!(mps2_pp, alg=Orthogonalize(trunc=trunc))
						end
					end
				end
				# println(norm(mps_pp), " ", norm(mps2_pp), " ", dot(mps_pp, mps2_pp), " ", distance(mps_pp, mps2_pp))
				@test distance(mps_pp, mps2_pp) / norm(mps2_pp) < tol

				mps2_pm = nothing
				for i in 1:lattice.k
					for j in 1:lattice.k
						pos1, pos2 = index(lattice, i, branch=:+, conj=true, band=band), index(lattice, j, branch=:-, conj=false, band=band)
						if corr.G₊₋[i, j] != 0
							t = GTerm(pos1, pos2, coeff=corr.G₊₋[i, j])
							if !isnothing(mps2_pm)
								mps2_pm += t * vc
							else
								mps2_pm = t * vc
							end
							canonicalize!(mps2_pm, alg=Orthogonalize(trunc=trunc))		
						end
					end
				end
				# println(norm(mps_pm), " ", norm(mps2_pm), " ",  dot(mps_pm, mps2_pm), " ", distance(mps_pm, mps2_pm))
				@test distance(mps_pm, mps2_pm) / norm(mps2_pm) < tol

				mps2_mp = nothing
				for i in 1:lattice.k
					for j in 1:lattice.k
						pos1, pos2 = index(lattice, i, branch=:-, conj=true, band=band), index(lattice, j, branch=:+, conj=false, band=band)
						if corr.G₋₊[i, j] != 0
							t = GTerm(pos1, pos2, coeff=corr.G₋₊[i, j])
							if !isnothing(mps2_mp)
								mps2_mp += t * vc
							else
								mps2_mp = t * vc
							end
							canonicalize!(mps2_mp, alg=Orthogonalize(trunc=trunc))
						end
					end
				end
				# println(norm(mps_mp), " ", norm(mps2_mp), " ", dot(mps_mp, mps2_mp), " ", distance(mps_mp, mps2_mp))
				@test distance(mps_mp, mps2_mp) / norm(mps2_mp) < tol				

				mps2_mm = nothing
				for i in 1:lattice.k
					for j in 1:lattice.k
						pos1, pos2 = index(lattice, i, branch=:-, conj=true, band=band), index(lattice, j, branch=:-, conj=false, band=band)
						if corr.G₋₋[i, j] != 0
							t = GTerm(pos1, pos2, coeff=corr.G₋₋[i, j])
							if !isnothing(mps2_mm)
								mps2_mm += t * vc
							else
								mps2_mm = t * vc
							end
							canonicalize!(mps2_mm, alg=Orthogonalize(trunc=trunc))
						end
					end
				end
				# println(norm(mps_mm), " ", norm(mps2_mm), " ", dot(mps_mm, mps2_mm), " ", distance(mps_mm, mps2_mm))
				@test distance(mps_mm, mps2_mm) / norm(mps2_mm) < tol

				dt = 0.01

				mps_pp = dt * mps_pp + vacuumstate(lattice)
				mps_pm = dt * mps_pm + vacuumstate(lattice)
				mps_mp = dt * mps_mp + vacuumstate(lattice)
				mps_mm = dt * mps_mm + vacuumstate(lattice)


				h_pp, h_pm, h_mp, h_mm = influenceoperatorexponential(lattice, corr, dt, WII(), band=band, algexpan=algexpan)
				mps2_pp = h_pp * vacuumstate(lattice)
				mps2_pm = h_pm * vacuumstate(lattice)
				mps2_mp = h_mp * vacuumstate(lattice)
				mps2_mm = h_mm * vacuumstate(lattice)

				@test distance(mps_pp, mps2_pp) / norm(mps2_pp) < dt
				@test distance(mps_pm, mps2_pm) / norm(mps2_pm) < dt
				@test distance(mps_mp, mps2_mp) / norm(mps2_mp) < dt
				@test distance(mps_mm, mps2_mm) / norm(mps2_mm) < dt

				mps0 = mult!(mps2_pp, mps2_pm, trunc=trunc)
				mps0 = mult!(mps0, mps2_mp, trunc=trunc)
				mps0 = mult!(mps0, mps2_mm, trunc=trunc)

				for algmult in (SVDCompression(D=50), DMRG1(D=50), DMRG2(trunc=truncdimcutoff(D=50,ϵ=1.0e-6)))
					mps1 = differentialinfluencefunctional(lattice, corr, dt, WII(), algmult, band=band, algexpan=algexpan)
					_n = norm(mps1)
					mps2 = differentialinfluencefunctional(lattice, corr, dt, WI(), algmult, band=band, algexpan=algexpan)
					mps3 = differentialinfluencefunctional(lattice, corr, dt, ComplexStepper(WI()), algmult, band=band, algexpan=algexpan)
					mps4 = differentialinfluencefunctional(lattice, corr, dt, ComplexStepper(WII()), algmult, band=band, algexpan=algexpan)
					@test distance(mps1, mps0) / _n < dt
					@test distance(mps2, mps0) / _n < dt
					@test distance(mps3, mps0) / _n < dt
					@test distance(mps4, mps0) / _n < dt
				end
			end

		end
	end

end

