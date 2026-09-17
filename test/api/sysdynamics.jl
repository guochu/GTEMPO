@testset "API: sysdynamics" begin
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	model = AndersonIM(U=1, ϵ_d=-0.5)

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
	ρ = fock_thermalstate(model, β)
	@test tr(ρ.data) ≈ 1 atol=1.0e-10
end

@testset "API: impurity hamiltonians & steppers" begin
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)
	model = AndersonIM(U=1, ϵ_d=-0.5)

	# --- predefined models: IRLM / KanamoriIM and the type hierarchy ---
	irlm = IRLM(μ=-0.5, J=0.3, U=1.0)
	@test irlm isa ImpurityHamiltonian
	@test irlm isa ConstImpurityHamiltonian && irlm isa AbstractImpurityHamiltonian
	@test num_bands(irlm) == 3
	kani = KanamoriIM(U=2.0, J=0.5, norb=1)
	@test kani isa ImpurityHamiltonian && num_bands(kani) == 2
	lat3 = GrassmannLattice(N=2, δτ=0.1, contour=:imag, bands=3)
	K3 = sysdynamics(lat3, irlm, trunc=trunc)
	@test length(K3) == length(lat3) && norm(K3) > 0

	# --- quenched and time-dependent hamiltonians (single band) ---
	hq = QuenchedImpurityHamiltonian([tunneling(1, 1, coeff=-0.5)], [tunneling(1, 1, coeff=0.5)]; bands=1)
	@test hq isa ConstImpurityHamiltonian
	htt = [TdImpurityOp([tunneling(1, 1)], t -> 0.1 * cos(t); bands=1)]
	hTd = TdImpurityHamiltonian([tunneling(1, 1, coeff=-0.5)], [tunneling(1, 1, coeff=0.5)], htt; bands=1)
	@test hTd isa AbstractTdImpurityHamiltonian
	@test hTd(0.7) isa ImpurityHamiltonian          # time-dependent part evaluated at t
	lat_q = GrassmannLattice(N=3, δt=0.05, contour=:real)
	Kq = sysdynamics(lat_q, hq, trunc=trunc)
	@test length(Kq) == length(lat_q) && norm(Kq) > 0
	lat1 = GrassmannLattice(N=3, δτ=0.1, contour=:imag)
	Ktd = sysdynamics(lat1, hTd, trunc=trunc)
	@test length(Ktd) == length(lat1) && norm(Ktd) > 0
	# the _fast variants require a constant (ConstImpurityHamiltonian) model
	@test_throws MethodError sysdynamics_fast(lat1, hTd, trunc=trunc)

	# --- branch-wise construction: forward (+ backward + imaginary) == full ---
	# (the analytic AndersonIM implementations with U ≠ 0 need a two-band lattice)
	lat = GrassmannLattice(N=3, δt=0.05, contour=:real, bands=2)
	Kp = sysdynamics(lat, model, branch=:+, trunc=trunc)
	@test distance(Kp, sysdynamics_forward!(vacuumstate(lat), lat, model, trunc=trunc)) < 1.0e-10
	Kfb = sysdynamics_forward!(vacuumstate(lat), lat, model, trunc=trunc)
	Kfb = sysdynamics_backward!(Kfb, lat, model, trunc=trunc)
	@test distance(Kfb, sysdynamics(lat, model, trunc=trunc)) < 1.0e-10
	# bare dynamics + bulk connection reproduces the full propagator
	bK = baresysdynamics_forward!(vacuumstate(lat), lat, model, trunc=trunc)
	bK = baresysdynamics_backward!(bK, lat, model, trunc=trunc)
	for band in 1:lat.bands
		bK = bulkconnection!(bK, lat, band=band)
	end
	@test distance(bK, sysdynamics(lat, model, trunc=trunc)) < 1.0e-8
	@test distance(baresysdynamics_fast(lat, model, trunc=trunc), baresysdynamics(lat, model, trunc=trunc)) < 1.0e-6
	# imaginary time
	lati = GrassmannLattice(N=3, δτ=0.1, contour=:imag, bands=2)
	Ki = sysdynamics_imaginary!(vacuumstate(lati), lati, model, trunc=trunc)
	@test distance(Ki, sysdynamics(lati, model, trunc=trunc)) < 1.0e-10
	bKi = baresysdynamics_imaginary!(vacuumstate(lati), lati, model, trunc=trunc)
	for band in 1:lati.bands
		bKi = bulkconnection!(bKi, lati, band=band)
	end
	@test distance(bKi, sysdynamics(lati, model, trunc=trunc)) / norm(Ki) < 1.0e-6
	# mixed time
	latm = GrassmannLattice(Nt=2, δt=0.05, Nτ=3, δτ=0.1, contour=:mixed, bands=2)
	Km = sysdynamics_forward!(vacuumstate(latm), latm, model, trunc=trunc)
	Km = sysdynamics_backward!(Km, latm, model, trunc=trunc)
	Km = sysdynamics_imaginary!(Km, latm, model, trunc=trunc)
	@test distance(Km, sysdynamics(latm, model, trunc=trunc)) < 1.0e-10

	# --- impurity initial states (real-time lattice, Fock-space operators) ---
	ρ0 = FockMatrix([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0])
	Ks = sysinitialstate(lat, ρ0; trunc=trunc)
	@test Ks isa GrassmannMPS && length(Ks) == length(lat)
	Ks2 = sysinitialstate!(vacuumstate(lat), lat, ρ0; trunc=trunc)
	@test distance(Ks, Ks2) < 1.0e-12
	Kt = systhermalstate(lat, model; trunc=trunc, β=1.0)
	@test Kt isa GrassmannMPS && length(Kt) == length(lat)
	Kt2 = systhermalstate!(vacuumstate(lat), lat, model, trunc=trunc, β=1.0)
	@test distance(Kt, Kt2) < 1.0e-12

	# --- stepwise (online) evolution: makestep + sysdynamicsstepper! ---
	# the stepped propagator agrees with the static construction up to the
	# coherent-state overlap factors carried in the overall normalization
	# (they cancel in normalized observables; see the validated stepwise
	# recipe in test/normalbath/continuousbath/stepwise.jl)
	lat_o = GrassmannLattice(N=3, δt=0.05, contour=:real, bands=2)
	lattice = similar(lat_o, N=0)
	mpsK = vacuumstate(lattice)
	while lattice.N < lat_o.N
		lattice, mpsK = makestep(lattice, mpsK)
		mpsK = sysdynamicsstepper!(mpsK, lattice, model, trunc=trunc)
	end
	@test length(mpsK) == length(lat_o) && norm(mpsK) > 0
	Kref = sysdynamics(lat_o, model, trunc=trunc)
	mpsK = mpsK * (norm(Kref) / norm(mpsK))
	@test distance(mpsK, Kref) / norm(Kref) < 1.0e-6

	# --- deprecated (trotterized) dynamics ---
	# (legacy convention: internally consistent, but its propagator differs
	# from the exact fock_propagator construction, so only family-internal
	# equalities and the bulk-connection closure are asserted)
	h1 = ImpurityHamiltonian(bands=1)
	push!(h1, tunneling(1, 1, coeff=-0.5))
	lat1i = GrassmannLattice(N=3, δτ=0.1, contour=:imag)
	Kdep = sysdynamics_deprecated(lat1i, h1, trunc=trunc)
	Kdep2 = sysdynamics_imaginary_deprecated!(vacuumstate(lat1i), lat1i, h1, trunc=trunc)
	@test distance(Kdep, Kdep2) < 1.0e-10
	Kdep3 = sysdynamics_deprecated!(vacuumstate(lat1i), lat1i, h1, trunc=trunc)
	@test distance(Kdep, Kdep3) < 1.0e-10
	bd1 = baresysdynamics_deprecated(lat1i, h1, trunc=trunc)
	bd2 = baresysdynamics_deprecated!(vacuumstate(lat1i), lat1i, h1, trunc=trunc)
	@test distance(bd1, bd2) < 1.0e-10
	bd1 = bulkconnection!(bd1, lat1i, band=1)
	@test distance(bd1, Kdep) / norm(Kdep) < 1.0e-8
	# real-time deprecated branches
	lat1r = GrassmannLattice(N=3, δt=0.05, contour=:real)
	Kdr = sysdynamics_deprecated(lat1r, h1, trunc=trunc)
	Kdr2 = sysdynamics_forward_deprecated!(vacuumstate(lat1r), lat1r, h1, trunc=trunc)
	Kdr2 = sysdynamics_backward_deprecated!(Kdr2, lat1r, h1, trunc=trunc)
	@test norm(Kdr) > 0
	@test distance(Kdr, Kdr2) < 1.0e-10
end
