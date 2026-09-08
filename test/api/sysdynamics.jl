@testset "API: sysdynamics" begin
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	model = AndersonIM(U=1, μ=-0.5)

	# imaginary time
	lat = GrassmannLattice(N=4, δτ=0.1, contour=:imag, bands=2)
	K = sysdynamics(lat, model, trunc=trunc)
	@test K isa GrassmannMPS && length(K) == length(lat) && norm(K) > 0
	for band in 1:lat.bands
		K = boundarycondition!(K, lat, band=band, trunc=trunc)
	end
	# baresysdynamics + bulkconnection reproduces sysdynamics
	K0 = baresysdynamics(lat, model, trunc=trunc)
	for band in 1:2
		K0 = bulkconnection!(K0, lat, band=band)
	end
	@test distance(K0, sysdynamics(lat, model, trunc=trunc)) / norm(K) < 1.0e-8

	# real time
	lat = GrassmannLattice(N=4, δt=0.05, contour=:real, bands=2)
	K = sysdynamics(lat, model, trunc=trunc)
	@test length(K) == length(lat) && norm(K) > 0
	for band in 1:lat.bands
		K = boundarycondition!(K, lat, band=band, trunc=trunc)
	end
	K = systhermalstate!(K, lat, model, trunc=trunc, β=1.0)
	# single-branch construction
	Kp = sysdynamics(lat, model, branch=:+, trunc=trunc)
	@test length(Kp) == length(lat)
	# fast variant is exact up to the truncation used by sysdynamics
	@test distance(sysdynamics_fast(lat, model, trunc=trunc), sysdynamics(lat, model, trunc=trunc)) < 1.0e-4

	# mixed time
	lat = GrassmannLattice(Nt=3, δt=0.05, Nτ=4, δτ=0.1, contour=:mixed, bands=2)
	K = sysdynamics(lat, model, trunc=trunc)
	@test length(K) == length(lat) && norm(K) > 0
	for band in 1:lat.bands
		K = boundarycondition!(K, lat, band=band, trunc=trunc)
	end
	Kτ = sysdynamics(lat, model, branch=:τ, trunc=trunc)
	@test length(Kτ) == length(lat)

	# generic ImpurityHamiltonian gives the same dynamics as AndersonIM
	h = ImpurityHamiltonian(bands=2)
	push!(h, interaction(1, 2, 2, 1, coeff=1))
	push!(h, tunneling(1, 1, coeff=-0.5))
	push!(h, tunneling(2, 2, coeff=-0.5))
	lat = GrassmannLattice(N=4, δτ=0.1, contour=:imag, bands=2)
	@test distance(sysdynamics(lat, h, trunc=trunc), sysdynamics(lat, model, trunc=trunc)) < 1.0e-5

	# thermal state on the impurity: partition function over the 2-band space
	β = 1.7
	ρ = fock_thermalstate(model, β, 2)
	@test tr(ρ.data) ≈ 1 atol=1.0e-10
end
