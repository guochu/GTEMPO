using LinearAlgebra: normalize!

@testset "API: normalize!" begin
	psi = randomgmps(Float64, 8, D=4)
	setscaling!(psi, 3.7)
	normalize!(psi)
	@test scaling(psi) ≈ 1.0 atol = 1.0e-12
	@test norm(psi[1]) ≈ 1.0 atol = 1.0e-10
end

@testset "API: integrateband / integratebands" begin
	lat3 = GrassmannLattice(N=2, δτ=0.1, bands=3, contour=:imag)
	lat2 = similar(lat3, bands=2)
	lat1 = similar(lat3, bands=1)
	x = randomgmps(Float64, length(lat3), D=4)
	Z = integrate(lat3, x)

	# integrating out one band and then the rest reproduces the full integral
	for b in 1:3
		xb = integrateband(lat3, x; band=b)
		@test length(xb) == length(lat2)
		@test relerr(Z, integrate(lat2, xb)) < 1.0e-8
	end

	# integratebands over two bands equals the manual cascade (descending bands)
	x13 = integratebands(lat3, x, (1, 3))
	@test length(x13) == length(lat1)
	@test relerr(Z, integrate(lat1, x13)) < 1.0e-8
	xm = integrateband(similar(lat3, bands=2), integrateband(lat3, x; band=3); band=1)
	@test relerr(Z, integrate(lat1, xm)) < 1.0e-8
	@test relerr(x13, xm) < 1.0e-8
end

@testset "API: multintegrateband" begin
	lat3 = GrassmannLattice(N=2, δτ=0.1, bands=3, contour=:imag)
	lat2 = similar(lat3, bands=2)
	trunc = truncdimcutoff(D=64, ϵ=1.0e-10)
	x = randomgmps(Float64, length(lat3), D=4)
	y = randomgmps(Float64, length(lat3), D=4)
	Z = integrate(lat3, x, y)

	xy1 = multintegrateband(lat3, x, y; trunc=trunc)
	xy2 = multintegrateband(lat3, x, y, SVDCompression(trunc))
	@test length(xy1) == length(lat2) && length(xy2) == length(lat2)
	@test relerr(Z, integrate(lat2, xy1)) < 1.0e-8
	@test relerr(Z, integrate(lat2, xy2)) < 1.0e-8
end

@testset "API: partialintegrate" begin
	# partialintegrate contracts several GMPS at once (2 or more are required)
	lat3 = GrassmannLattice(N=2, δτ=0.1, bands=3, contour=:imag)
	lat2 = similar(lat3, bands=2)
	lat1 = similar(lat3, bands=1)
	trunc = truncdimcutoff(D=64, ϵ=1.0e-10)
	x = randomgmps(Float64, length(lat3), D=4)
	y = randomgmps(Float64, length(lat3), D=4)
	z = randomgmps(Float64, length(lat3), D=4)
	Zxy = integrate(lat3, x, y)

	# both the SVD compression and the iterative DMRG1 algorithm are supported
	for alg in (SVDCompression(trunc), DMRG1(trunc))
		# two GMPS, integrating out a single band
		xyp = partialintegrate(lat3, x, y; branchs=(:τ,), bands=(1,), alg=alg)
		@test length(xyp) == length(lat2)
		@test relerr(Zxy, integrate(lat2, xyp)) < 1.0e-8

		# two GMPS, integrating out two bands at once
		xyp13 = partialintegrate(lat3, x, y; branchs=(:τ,), bands=(1, 3), alg=alg)
		@test length(xyp13) == length(lat1)
		@test relerr(Zxy, integrate(lat1, xyp13)) < 1.0e-8

		# three GMPS
		xyzp = partialintegrate(lat3, x, y, z; branchs=(:τ,), bands=(1,), alg=alg)
		@test length(xyzp) == length(lat2)
		@test relerr(integrate(lat3, x, y, z), integrate(lat2, xyzp)) < 1.0e-8
	end
end

@testset "API: zipup integration, integrate(A,B,...) = integrate(A*B...)" begin
	# the naive product of k random GMPS with D=2 stays small enough to be exact,
	# while the zipup algorithms contract them step by step
	for (contour, lkw, T) in ((:imag, (N=2, δτ=0.1), Float64), (:real, (N=2, δt=0.05), ComplexF64))
		lat = GrassmannLattice(; contour=contour, lkw...)
		L = length(lat)
		xs = [randomgmps(T, L, D=2) for _ in 1:5]
		for k in (2, 5)
			z_prod = integrate(lat, *(xs[1:k]...))
			z_zip = integrate(lat, xs[1:k]...)
			@test relerr(z_prod, z_zip) < 1.0e-8
		end
		# the boundary-MPS variant agrees with the exact zipup
		trunc = truncdimcutoff(D=64, ϵ=1.0e-10)
		z_exact = integrate(lat, xs[1:3]...)
		z_bmps = integrate(lat, xs[1:3]...; alg=BMPSIntegrate(trunc))
		@test relerr(z_exact, z_bmps) < 1.0e-8
	end
end
