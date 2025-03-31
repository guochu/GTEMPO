println("------------------------------------")
println("|             FockLattice          |")
println("------------------------------------")


@testset "FockLattice: imaginary time M1N1" begin
	lattice = FockLattice(N=1, δτ=0.1, contour=:imag, ordering=M1N1())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, ImagFockLattice)
	@test lattice.ordering == M1N1()
	@test scalartype(lattice) == Float64
	@test length(lattice) == 1
	@test lattice.N == 1
	@test lattice.β == 0.1
	@test lattice.bands == 1
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.1
	@test lattice.T == 10
	@test index(lattice, 1) == 1

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 1 atol = 1.0e-6

	lattice = FockLattice(N=2, δτ=0.05, bands=2, contour=:imag, ordering=M1N1())
	@test isa(lattice, ImagFockLattice)
	@test length(lattice) == 4
	@test lattice.N == 2
	@test lattice.β == 0.1
	@test lattice.bands == 2
	@test lattice.δτ == 0.05
	@test lattice.τs == 0:0.05:0.1
	@test lattice.T == 10

	@test index(lattice, 2, band=1) == 1
	@test index(lattice, 2, band=2) == 2
	@test index(lattice, 1, band=1) == 3
	@test index(lattice, 1, band=2) == 4

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 1 atol = 1.0e-6

	lattice = FockLattice(N=2, δτ=0.1, bands=3, contour=:imag, ordering=M1N1())
	@test length(lattice) == 6
	@test lattice.N == 2
	@test lattice.β == 0.2
	@test lattice.bands == 3
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.2	

	@test index(lattice, 2, band=1) == 1
	@test index(lattice, 2, band=2) == 2
	@test index(lattice, 2, band=3) == 3
	@test index(lattice, 1, band=1) == 4
	@test index(lattice, 1, band=2) == 5
	@test index(lattice, 1, band=3) == 6

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 1 atol = 1.0e-6
end


@testset "FockLattice: real time M1m1N1n1" begin
	lattice = FockLattice(N=1, δt=0.1, contour=:real, ordering=M1m1N1n1())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, RealFockLattice)
	@test lattice.ordering == M1m1N1n1()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 2
	@test lattice.N == 1
	@test lattice.t == 0.1
	@test lattice.bands == 1
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 1, branch=:+) == 1
	@test index(lattice, 1, branch=:-) == 2

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 1 atol = 1.0e-6

	lattice = FockLattice(N=2, δt=0.05, bands=2, contour=:real, ordering=M1m1N1n1())
	@test isa(lattice, RealFockLattice)
	@test length(lattice) == 8
	@test lattice.N == 2
	@test lattice.t == 0.1
	@test lattice.bands == 2
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1

	@test index(lattice, 2, band=1, branch=:+) == 1
	@test index(lattice, 2, band=1, branch=:-) == 2
	@test index(lattice, 2, band=2, branch=:+) == 3
	@test index(lattice, 2, band=2, branch=:-) == 4	

	@test index(lattice, 1, band=1, branch=:+) == 5
	@test index(lattice, 1, band=1, branch=:-) == 6
	@test index(lattice, 1, band=2, branch=:+) == 7
	@test index(lattice, 1, band=2, branch=:-) == 8	

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 1 atol = 1.0e-6

	lattice = FockLattice(N=2, δt=0.1, bands=3, contour=:real, ordering=M1m1N1n1())
	@test length(lattice) == 12
	@test lattice.N == 2
	@test lattice.t == 0.2
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.2	

	@test index(lattice, 2, band=1, branch=:+) == 1
	@test index(lattice, 2, band=1, branch=:-) == 2
	@test index(lattice, 2, band=2, branch=:+) == 3
	@test index(lattice, 2, band=2, branch=:-) == 4	
	@test index(lattice, 2, band=3, branch=:+) == 5
	@test index(lattice, 2, band=3, branch=:-) == 6	

	@test index(lattice, 1, band=1, branch=:+) == 7
	@test index(lattice, 1, band=1, branch=:-) == 8
	@test index(lattice, 1, band=2, branch=:+) == 9
	@test index(lattice, 1, band=2, branch=:-) == 10	
	@test index(lattice, 1, band=3, branch=:+) == 11
	@test index(lattice, 1, band=3, branch=:-) == 12

	mps = vacuumstate(lattice)
	@test integrate(mps) ≈ 1 atol = 1.0e-6
end

@testset "FockLattice: mixed time M1N1_m1M1n1N1m2M2n2N2" begin
	# one band
	lattice = FockLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=M1N1_m1M1n1N1m2M2n2N2())
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, MixedFockLattice)
	@test lattice.ordering == M1N1_m1M1n1N1m2M2n2N2()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 6
	@test lattice.Nt == 2
	@test lattice.Nτ == 2
	@test lattice.t == 0.1
	@test lattice.β == 0.2
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test lattice.τs == 0:0.1:0.2

	# imaginary time axis
	@test index(lattice, 2, branch=:τ) == 1
	@test index(lattice, 1, branch=:τ) == 2

	# real time axis
	@test index(lattice, 1, branch=:-) == 3
	@test index(lattice, 1, branch=:+) == 4

	@test index(lattice, 2, branch=:-) == 5
	@test index(lattice, 2, branch=:+) == 6

	# two bands
	lattice = FockLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, contour=:mixed, bands=2, ordering=M1N1_m1M1n1N1m2M2n2N2())
	@test isa(lattice, MixedFockLattice)
	@test length(lattice) == 8
	@test lattice.Nt == 1
	@test lattice.Nτ == 2
	@test lattice.t == 0.1
	@test lattice.β == 0.2
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test lattice.τs == 0:0.1:0.2


	# imaginary time axis
	@test index(lattice, 2, branch=:τ, band=1) == 1
	@test index(lattice, 2, branch=:τ, band=2) == 2
	@test index(lattice, 1, branch=:τ, band=1) == 3
	@test index(lattice, 1, branch=:τ, band=2) == 4

	# real time axis
	@test index(lattice, 1, branch=:-, band=1) == 5
	@test index(lattice, 1, branch=:+, band=1) == 6
	@test index(lattice, 1, branch=:-, band=2) == 7
	@test index(lattice, 1, branch=:+, band=2) == 8

	# three bands
	lattice = FockLattice(Nt=1, δt=0.1, Nτ=1, δτ=0.1, contour=:mixed, bands=3, ordering=M1N1_m1M1n1N1m2M2n2N2())
	@test isa(lattice, MixedFockLattice)
	@test length(lattice) == 9
	@test lattice.Nt == 1
	@test lattice.Nτ == 1
	@test lattice.t == 0.1
	@test lattice.β == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test lattice.τs == 0:0.1:0.1


	# imaginary time axis
	@test index(lattice, 1, branch=:τ, band=1) == 1
	@test index(lattice, 1, branch=:τ, band=2) == 2
	@test index(lattice, 1, branch=:τ, band=3) == 3

	# real time axis
	@test index(lattice, 1, branch=:-, band=1) == 4
	@test index(lattice, 1, branch=:+, band=1) == 5
	@test index(lattice, 1, branch=:-, band=2) == 6
	@test index(lattice, 1, branch=:+, band=2) == 7
	@test index(lattice, 1, branch=:-, band=3) == 8
	@test index(lattice, 1, branch=:+, band=3) == 9
end

