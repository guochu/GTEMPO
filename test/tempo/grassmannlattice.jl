println("------------------------------------")
println("|          GrassmannLattice        |")
println("------------------------------------")

@testset "GrassmannLattice: imaginary time A1B1B1A1" begin
	# a\bar{a} a_2\bar{a}_2 a_1\bar{a}_1
	lattice = GrassmannLattice(N=1, δτ=0.1, contour=:imag, ordering=ABBA())
	@test isa(ConjugationStyle(lattice), GeneralConjugation)
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, ImagGrassmannLattice)
	@test lattice.ordering == ABBA()
	@test scalartype(lattice) == Float64
	@test length(lattice) == 6
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.β == 0.1
	@test lattice.bands == 1
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.1
	@test lattice.T == 10
	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2
	@test index(lattice, 2, conj=false) == 3
	@test index(lattice, 2, conj=true) == 4
	@test index(lattice, 1, conj=false) == 5
	@test index(lattice, 1, conj=true) == 6

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# ab\bar{b}\bar{a} a_3b_3\bar{b}_3\bar{a}_3 a_2b_2\bar{b}_2\bar{a}_2 a_1b_1\bar{b}_1\bar{a}_1
	lattice = GrassmannLattice(N=2, δτ=0.05, bands=2, contour=:imag, ordering=ABBA())
	@test isa(lattice, ImagGrassmannLattice)
	@test length(lattice) == 16
	@test lattice.N == 2
	@test lattice.k == 3
	@test lattice.β == 0.1
	@test lattice.bands == 2
	@test lattice.δτ == 0.05
	@test lattice.τs == 0:0.05:0.1
	@test lattice.T == 10

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 2
	@test index(lattice, 0, conj=true, band=1) == 4
	@test index(lattice, 0, conj=true, band=2) == 3

	@test index(lattice, 3, conj=false, band=1) == 5
	@test index(lattice, 3, conj=false, band=2) == 6
	@test index(lattice, 3, conj=true, band=1) == 8
	@test index(lattice, 3, conj=true, band=2) == 7

	@test index(lattice, 2, conj=false, band=1) == 9
	@test index(lattice, 2, conj=false, band=2) == 10
	@test index(lattice, 2, conj=true, band=1) == 12
	@test index(lattice, 2, conj=true, band=2) == 11

	@test index(lattice, 1, conj=false, band=1) == 13
	@test index(lattice, 1, conj=false, band=2) == 14
	@test index(lattice, 1, conj=true, band=1) == 16
	@test index(lattice, 1, conj=true, band=2) == 15

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# abc\bar{c}\bar{b}\bar{a} a_2b_2c_2\bar{c}_2\bar{b}_2\bar{a}_2 a_1b_1c_1\bar{c}_1\bar{b}_1\bar{a}_1
	lattice = GrassmannLattice(N=1, δτ=0.1, bands=3, contour=:imag, ordering=ABBA())
	@test length(lattice) == 18
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.β == 0.1
	@test lattice.bands == 3
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.1	
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 2
	@test index(lattice, 0, conj=false, band=3) == 3
	@test index(lattice, 0, conj=true, band=1) == 6
	@test index(lattice, 0, conj=true, band=2) == 5
	@test index(lattice, 0, conj=true, band=3) == 4
	@test index(lattice, 2, conj=false, band=1) == 7
	@test index(lattice, 2, conj=false, band=2) == 8
	@test index(lattice, 2, conj=false, band=3) == 9
	@test index(lattice, 2, conj=true, band=1) == 12
	@test index(lattice, 2, conj=true, band=2) == 11
	@test index(lattice, 2, conj=true, band=3) == 10
	@test index(lattice, 1, conj=false, band=1) == 13
	@test index(lattice, 1, conj=false, band=2) == 14
	@test index(lattice, 1, conj=false, band=3) == 15
	@test index(lattice, 1, conj=true, band=1) == 18
	@test index(lattice, 1, conj=true, band=2) == 17
	@test index(lattice, 1, conj=true, band=3) == 16

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

@testset "GrassmannLattice: imaginary time A1A1B1B1" begin
	# a\bar{a} a_2\bar{a}_2 a_1\bar{a}_1
	lattice = GrassmannLattice(N=1, δτ=0.1, contour=:imag, ordering=AABB())
	@test isa(ConjugationStyle(lattice), AdjacentConjugation)
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, ImagGrassmannLattice)
	@test lattice.ordering == AABB()
	@test scalartype(lattice) == Float64
	@test length(lattice) == 6
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.β == 0.1
	@test lattice.bands == 1
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.1
	@test lattice.T == 10
	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2
	@test index(lattice, 2, conj=false) == 3
	@test index(lattice, 2, conj=true) == 4
	@test index(lattice, 1, conj=false) == 5
	@test index(lattice, 1, conj=true) == 6

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# a\bar{a}b\bar{b} a_3\bar{a}_3b_3\bar{b}_3 a_2\bar{a}_2b_2\bar{b}_2 a_1\bar{a}_1b_1\bar{b}_1
	lattice = GrassmannLattice(N=2, δτ=0.05, bands=2, contour=:imag, ordering=AABB())
	@test isa(lattice, ImagGrassmannLattice)
	@test length(lattice) == 16
	@test lattice.N == 2
	@test lattice.k == 3
	@test lattice.β == 0.1
	@test lattice.bands == 2
	@test lattice.δτ == 0.05
	@test lattice.τs == 0:0.05:0.1
	@test lattice.T == 10

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=true, band=2) == 4

	@test index(lattice, 3, conj=false, band=1) == 5
	@test index(lattice, 3, conj=false, band=2) == 7
	@test index(lattice, 3, conj=true, band=1) == 6
	@test index(lattice, 3, conj=true, band=2) == 8

	@test index(lattice, 2, conj=false, band=1) == 9
	@test index(lattice, 2, conj=false, band=2) == 11
	@test index(lattice, 2, conj=true, band=1) == 10
	@test index(lattice, 2, conj=true, band=2) == 12

	@test index(lattice, 1, conj=false, band=1) == 13
	@test index(lattice, 1, conj=false, band=2) == 15
	@test index(lattice, 1, conj=true, band=1) == 14
	@test index(lattice, 1, conj=true, band=2) == 16

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# abc\bar{c}\bar{b}\bar{a} a_2b_2c_2\bar{c}_2\bar{b}_2\bar{a}_2 a_1b_1c_1\bar{c}_1\bar{b}_1\bar{a}_1
	lattice = GrassmannLattice(N=1, δτ=0.1, bands=3, contour=:imag, ordering=AABB())
	@test length(lattice) == 18
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.β == 0.1
	@test lattice.bands == 3
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.1	
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4
	@test index(lattice, 0, conj=false, band=3) == 5
	@test index(lattice, 0, conj=true, band=3) == 6

	@test index(lattice, 2, conj=false, band=1) == 7
	@test index(lattice, 2, conj=true, band=1) == 8
	@test index(lattice, 2, conj=false, band=2) == 9
	@test index(lattice, 2, conj=true, band=2) == 10
	@test index(lattice, 2, conj=false, band=3) == 11
	@test index(lattice, 2, conj=true, band=3) == 12
	@test index(lattice, 1, conj=false, band=1) == 13
	@test index(lattice, 1, conj=true, band=1) == 14
	@test index(lattice, 1, conj=false, band=2) == 15
	@test index(lattice, 1, conj=true, band=2) == 16
	@test index(lattice, 1, conj=false, band=3) == 17
	@test index(lattice, 1, conj=true, band=3) == 18

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

@testset "GrassmannLattice: imaginary time A2A2A1A1B2B2B1B1" begin
	# one band
	lattice = GrassmannLattice(N=1, δτ=0.1, contour=:imag, ordering=A2A2A1A1B2B2B1B1())
	@test isa(ConjugationStyle(lattice), AdjacentConjugation)
	@test isa(LayoutStyle(lattice), BandLocalLayout)
	@test isa(lattice, ImagGrassmannLattice)
	@test lattice.ordering == A2A2A1A1B2B2B1B1()
	@test scalartype(lattice) == Float64
	@test length(lattice) == 6
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.β == 0.1
	@test lattice.bands == 1
	@test lattice.δτ == 0.1
	@test lattice.τs == 0:0.1:0.1
	@test lattice.T == 10
	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2
	@test index(lattice, 2, conj=false) == 3
	@test index(lattice, 2, conj=true) == 4
	@test index(lattice, 1, conj=false) == 5
	@test index(lattice, 1, conj=true) == 6

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# two bands
	lattice = GrassmannLattice(N=1, δτ=0.05, bands=2, contour=:imag, ordering=A2A2A1A1B2B2B1B1())
	@test length(lattice) == 12
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.β == 0.05
	@test lattice.bands == 2
	@test lattice.δτ == 0.05
	@test lattice.τs == 0:0.05:0.05
	@test lattice.T == 20

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=true, band=2) == 4

	@test index(lattice, 2, conj=false, band=1) == 5
	@test index(lattice, 2, conj=true, band=1) == 6
	@test index(lattice, 1, conj=false, band=1) == 7
	@test index(lattice, 1, conj=true, band=1) == 8

	@test index(lattice, 2, conj=false, band=2) == 9
	@test index(lattice, 2, conj=true, band=2) == 10
	@test index(lattice, 1, conj=false, band=2) == 11
	@test index(lattice, 1, conj=true, band=2) == 12

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6


	# three bands
	lattice = GrassmannLattice(N=1, δτ=0.05, bands=3, contour=:imag, ordering=A2A2A1A1B2B2B1B1())
	@test length(lattice) == 18
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.β == 0.05
	@test lattice.bands == 3
	@test lattice.δτ == 0.05
	@test lattice.τs == 0:0.05:0.05
	@test lattice.T == 20

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4
	@test index(lattice, 0, conj=false, band=3) == 5
	@test index(lattice, 0, conj=true, band=3) == 6

	@test index(lattice, 2, conj=false, band=1) == 7
	@test index(lattice, 2, conj=true, band=1) == 8
	@test index(lattice, 1, conj=false, band=1) == 9
	@test index(lattice, 1, conj=true, band=1) == 10

	@test index(lattice, 2, conj=false, band=2) == 11
	@test index(lattice, 2, conj=true, band=2) == 12
	@test index(lattice, 1, conj=false, band=2) == 13
	@test index(lattice, 1, conj=true, band=2) == 14

	@test index(lattice, 2, conj=false, band=3) == 15
	@test index(lattice, 2, conj=true, band=3) == 16
	@test index(lattice, 1, conj=false, band=3) == 17
	@test index(lattice, 1, conj=true, band=3) == 18

	mps = GrassmannMPS(length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

@testset "GrassmannLattice: real time A1a1B1b1b1B1a1A1" begin
	# a\bar{a} a_3^+a_3^-\bar{a}_3^-\bar{a}_3^+ a_2^+a_2^-\bar{a}_2^-\bar{a}_2^+ a_1^+a_1^-\bar{a}_1^-\bar{a}_1^+
	lattice = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A1a1B1b1b1B1a1A1())
	@test isa(ConjugationStyle(lattice), GeneralConjugation)
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, RealGrassmannLattice)
	@test lattice.ordering == A1a1B1b1b1B1a1A1()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 14
	@test lattice.N == 2
	@test lattice.k == 3
	@test lattice.t == 0.1
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2
	@test index(lattice, 3, conj=false, forward=true) == 3
	@test index(lattice, 3, conj=false, forward=false) == 4
	@test index(lattice, 3, conj=true, forward=true) == 6
	@test index(lattice, 3, conj=true, forward=false) == 5
	@test index(lattice, 2, conj=false, forward=true) == 7
	@test index(lattice, 2, conj=false, forward=false) == 8
	@test index(lattice, 2, conj=true, forward=true) == 10
	@test index(lattice, 2, conj=true, forward=false) == 9	
	@test index(lattice, 1, conj=false, forward=true) == 11
	@test index(lattice, 1, conj=false, forward=false) == 12
	@test index(lattice, 1, conj=true, forward=true) == 14
	@test index(lattice, 1, conj=true, forward=false) == 13	

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# ab\bar{b}\bar{a} a_2^+a_2^-b_2^+b_2^-\bar{b}_2^-\bar{b}_2^+\bar{a}_2^-\bar{a}_2^+ a_1^+a_1^-b_1^+b_1^-\bar{b}_1^-\bar{b}_1^+\bar{a}_1^-\bar{a}_1^+
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=2, ordering=A1a1B1b1b1B1a1A1())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 20
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 2
	@test index(lattice, 0, conj=true, band=1) == 4
	@test index(lattice, 0, conj=true, band=2) == 3

	@test index(lattice, 2, conj=false, forward=true, band=1) == 5
	@test index(lattice, 2, conj=false, forward=false, band=1) == 6
	@test index(lattice, 2, conj=false, forward=true, band=2) == 7
	@test index(lattice, 2, conj=false, forward=false, band=2) == 8
	@test index(lattice, 2, conj=true, forward=true, band=1) == 12
	@test index(lattice, 2, conj=true, forward=false, band=1) == 11
	@test index(lattice, 2, conj=true, forward=true, band=2) == 10
	@test index(lattice, 2, conj=true, forward=false, band=2) == 9

	@test index(lattice, 1, conj=false, forward=true, band=1) == 13
	@test index(lattice, 1, conj=false, forward=false, band=1) == 14
	@test index(lattice, 1, conj=false, forward=true, band=2) == 15
	@test index(lattice, 1, conj=false, forward=false, band=2) == 16
	@test index(lattice, 1, conj=true, forward=true, band=1) == 20
	@test index(lattice, 1, conj=true, forward=false, band=1) == 19
	@test index(lattice, 1, conj=true, forward=true, band=2) == 18
	@test index(lattice, 1, conj=true, forward=false, band=2) == 17

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# abc\bar{c}\bar{b}\bar{a} a_2^+a_2^-b_2^+b_2^-c_2^+c_2^-\bar{c}_2^-\bar{c}_2^+\bar{b}_2^-\bar{b}_2^+\bar{a}_2^-\bar{a}_2^+ 
	# a_1^+a_1^-b_1^+b_1^-c_1^+c_1^-\bar{c}_1^-\bar{c}_1^+\bar{b}_1^-\bar{b}_1^+\bar{a}_1^-\bar{a}_1^+
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=3, ordering=A1a1B1b1b1B1a1A1())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 30
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 2
	@test index(lattice, 0, conj=false, band=3) == 3
	@test index(lattice, 0, conj=true, band=1) == 6
	@test index(lattice, 0, conj=true, band=2) == 5	
	@test index(lattice, 0, conj=true, band=3) == 4

	@test index(lattice, 2, conj=false, forward=true, band=1) == 7
	@test index(lattice, 2, conj=false, forward=false, band=1) == 8
	@test index(lattice, 2, conj=false, forward=true, band=2) == 9
	@test index(lattice, 2, conj=false, forward=false, band=2) == 10
	@test index(lattice, 2, conj=false, forward=true, band=3) == 11
	@test index(lattice, 2, conj=false, forward=false, band=3) == 12
	@test index(lattice, 2, conj=true, forward=true, band=1) == 18
	@test index(lattice, 2, conj=true, forward=false, band=1) == 17
	@test index(lattice, 2, conj=true, forward=true, band=2) == 16
	@test index(lattice, 2, conj=true, forward=false, band=2) == 15
	@test index(lattice, 2, conj=true, forward=true, band=3) == 14
	@test index(lattice, 2, conj=true, forward=false, band=3) == 13

	@test index(lattice, 1, conj=false, forward=true, band=1) == 19
	@test index(lattice, 1, conj=false, forward=false, band=1) == 20
	@test index(lattice, 1, conj=false, forward=true, band=2) == 21
	@test index(lattice, 1, conj=false, forward=false, band=2) == 22
	@test index(lattice, 1, conj=false, forward=true, band=3) == 23
	@test index(lattice, 1, conj=false, forward=false, band=3) == 24
	@test index(lattice, 1, conj=true, forward=true, band=1) == 30
	@test index(lattice, 1, conj=true, forward=false, band=1) == 29
	@test index(lattice, 1, conj=true, forward=true, band=2) == 28
	@test index(lattice, 1, conj=true, forward=false, band=2) == 27
	@test index(lattice, 1, conj=true, forward=true, band=3) == 26
	@test index(lattice, 1, conj=true, forward=false, band=3) == 25

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end


@testset "GrassmannLattice: real time A1A1a1a1B1B1b1b1" begin
	# a\bar{a} a_3^+\bar{a}_3^+a_3^-\bar{a}_3^- a_2^+\bar{a}_2^+a_2^-\bar{a}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-
	lattice = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A1A1a1a1B1B1b1b1())
	@test isa(ConjugationStyle(lattice), AdjacentConjugation)
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, RealGrassmannLattice)
	@test lattice.ordering == A1A1a1a1B1B1b1b1()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 14
	@test lattice.N == 2
	@test lattice.k == 3
	@test lattice.t == 0.1
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2
	@test index(lattice, 3, conj=false, forward=true) == 3
	@test index(lattice, 3, conj=true, forward=true) == 4
	@test index(lattice, 3, conj=false, forward=false) == 5
	@test index(lattice, 3, conj=true, forward=false) == 6
	@test index(lattice, 2, conj=false, forward=true) == 7
	@test index(lattice, 2, conj=true, forward=true) == 8
	@test index(lattice, 2, conj=false, forward=false) == 9
	@test index(lattice, 2, conj=true, forward=false) == 10
	@test index(lattice, 1, conj=false, forward=true) == 11
	@test index(lattice, 1, conj=true, forward=true) == 12
	@test index(lattice, 1, conj=false, forward=false) == 13
	@test index(lattice, 1, conj=true, forward=false) == 14

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# a\bar{a}b\bar{b} a_2^+\bar{a}_2^+a_2^-\bar{a}_2^-b_2^+\bar{b}_2^+b_2^-\bar{b}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-b_1^+\bar{b}_1^+b_1^-\bar{b}_1^-
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=2, ordering=A1A1a1a1B1B1b1b1())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 20
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4

	@test index(lattice, 2, conj=false, forward=true, band=1) == 5
	@test index(lattice, 2, conj=true, forward=true, band=1) == 6
	@test index(lattice, 2, conj=false, forward=false, band=1) == 7
	@test index(lattice, 2, conj=true, forward=false, band=1) == 8
	@test index(lattice, 2, conj=false, forward=true, band=2) == 9
	@test index(lattice, 2, conj=true, forward=true, band=2) == 10
	@test index(lattice, 2, conj=false, forward=false, band=2) == 11
	@test index(lattice, 2, conj=true, forward=false, band=2) == 12

	@test index(lattice, 1, conj=false, forward=true, band=1) == 13
	@test index(lattice, 1, conj=true, forward=true, band=1) == 14
	@test index(lattice, 1, conj=false, forward=false, band=1) == 15
	@test index(lattice, 1, conj=true, forward=false, band=1) == 16
	@test index(lattice, 1, conj=false, forward=true, band=2) == 17
	@test index(lattice, 1, conj=true, forward=true, band=2) == 18
	@test index(lattice, 1, conj=false, forward=false, band=2) == 19
	@test index(lattice, 1, conj=true, forward=false, band=2) == 20

	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# a\bar{a}b\bar{b}c\bar{c} a_2^+\bar{a}_2^+a_2^-\bar{a}_2^-b_2^+\bar{b}_2^+b_2^-\bar{b}_2^-c_2^+\bar{c}_2^+c_2^-\bar{c}_2^- 
	# a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-b_1^+\bar{b}_1^+b_1^-\bar{b}_1^-c_1^+\bar{c}_1^+c_1^-\bar{c}_1^-
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=3, ordering=A1A1a1a1B1B1b1b1())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 30
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4
	@test index(lattice, 0, conj=false, band=3) == 5
	@test index(lattice, 0, conj=true, band=3) == 6

	@test index(lattice, 2, conj=false, forward=true, band=1) == 7
	@test index(lattice, 2, conj=true, forward=true, band=1) == 8
	@test index(lattice, 2, conj=false, forward=false, band=1) == 9
	@test index(lattice, 2, conj=true, forward=false, band=1) == 10
	@test index(lattice, 2, conj=false, forward=true, band=2) == 11
	@test index(lattice, 2, conj=true, forward=true, band=2) == 12
	@test index(lattice, 2, conj=false, forward=false, band=2) == 13
	@test index(lattice, 2, conj=true, forward=false, band=2) == 14
	@test index(lattice, 2, conj=false, forward=true, band=3) == 15
	@test index(lattice, 2, conj=true, forward=true, band=3) == 16
	@test index(lattice, 2, conj=false, forward=false, band=3) == 17
	@test index(lattice, 2, conj=true, forward=false, band=3) == 18

	@test index(lattice, 1, conj=false, forward=true, band=1) == 19
	@test index(lattice, 1, conj=true, forward=true, band=1) == 20
	@test index(lattice, 1, conj=false, forward=false, band=1) == 21
	@test index(lattice, 1, conj=true, forward=false, band=1) == 22
	@test index(lattice, 1, conj=false, forward=true, band=2) == 23
	@test index(lattice, 1, conj=true, forward=true, band=2) == 24
	@test index(lattice, 1, conj=false, forward=false, band=2) == 25
	@test index(lattice, 1, conj=true, forward=false, band=2) == 26
	@test index(lattice, 1, conj=false, forward=true, band=3) == 27
	@test index(lattice, 1, conj=true, forward=true, band=3) == 28
	@test index(lattice, 1, conj=false, forward=false, band=3) == 29
	@test index(lattice, 1, conj=true, forward=false, band=3) == 30

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end



@testset "GrassmannLattice: real time A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1" begin
	# one band
	lattice = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1())
	@test isa(ConjugationStyle(lattice), AdjacentConjugation)
	@test isa(LayoutStyle(lattice), BandLocalLayout)
	@test isa(lattice, RealGrassmannLattice)
	@test lattice.ordering == A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 14
	@test lattice.N == 2
	@test lattice.k == 3
	@test lattice.t == 0.1
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2
	@test index(lattice, 3, conj=false, forward=true) == 3
	@test index(lattice, 3, conj=true, forward=true) == 4
	@test index(lattice, 2, conj=false, forward=true) == 5
	@test index(lattice, 2, conj=true, forward=true) == 6
	@test index(lattice, 1, conj=false, forward=true) == 7
	@test index(lattice, 1, conj=true, forward=true) == 8
	@test index(lattice, 3, conj=false, forward=false) == 9
	@test index(lattice, 3, conj=true, forward=false) == 10
	@test index(lattice, 2, conj=false, forward=false) == 11
	@test index(lattice, 2, conj=true, forward=false) == 12
	@test index(lattice, 1, conj=false, forward=false) == 13
	@test index(lattice, 1, conj=true, forward=false) == 14

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# two bands
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=2, ordering=A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 20
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4

	@test index(lattice, 2, conj=false, forward=true, band=1) == 5
	@test index(lattice, 2, conj=true, forward=true, band=1) == 6
	@test index(lattice, 1, conj=false, forward=true, band=1) == 7
	@test index(lattice, 1, conj=true, forward=true, band=1) == 8
	@test index(lattice, 2, conj=false, forward=false, band=1) == 9
	@test index(lattice, 2, conj=true, forward=false, band=1) == 10
	@test index(lattice, 1, conj=false, forward=false, band=1) == 11
	@test index(lattice, 1, conj=true, forward=false, band=1) == 12

	@test index(lattice, 2, conj=false, forward=true, band=2) == 13
	@test index(lattice, 2, conj=true, forward=true, band=2) == 14
	@test index(lattice, 1, conj=false, forward=true, band=2) == 15
	@test index(lattice, 1, conj=true, forward=true, band=2) == 16
	@test index(lattice, 2, conj=false, forward=false, band=2) == 17
	@test index(lattice, 2, conj=true, forward=false, band=2) == 18
	@test index(lattice, 1, conj=false, forward=false, band=2) == 19
	@test index(lattice, 1, conj=true, forward=false, band=2) == 20

	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6


	# two bands
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=3, ordering=A2A2A1A1a2a2a1a1B2B2B1B1b2b2b1b1())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 30
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4
	@test index(lattice, 0, conj=false, band=3) == 5
	@test index(lattice, 0, conj=true, band=3) == 6

	@test index(lattice, 2, conj=false, forward=true, band=1) == 7
	@test index(lattice, 2, conj=true, forward=true, band=1) == 8
	@test index(lattice, 1, conj=false, forward=true, band=1) == 9
	@test index(lattice, 1, conj=true, forward=true, band=1) == 10
	@test index(lattice, 2, conj=false, forward=false, band=1) == 11
	@test index(lattice, 2, conj=true, forward=false, band=1) == 12
	@test index(lattice, 1, conj=false, forward=false, band=1) == 13
	@test index(lattice, 1, conj=true, forward=false, band=1) == 14

	@test index(lattice, 2, conj=false, forward=true, band=2) == 15
	@test index(lattice, 2, conj=true, forward=true, band=2) == 16
	@test index(lattice, 1, conj=false, forward=true, band=2) == 17
	@test index(lattice, 1, conj=true, forward=true, band=2) == 18
	@test index(lattice, 2, conj=false, forward=false, band=2) == 19
	@test index(lattice, 2, conj=true, forward=false, band=2) == 20
	@test index(lattice, 1, conj=false, forward=false, band=2) == 21
	@test index(lattice, 1, conj=true, forward=false, band=2) == 22

	@test index(lattice, 2, conj=false, forward=true, band=3) == 23
	@test index(lattice, 2, conj=true, forward=true, band=3) == 24
	@test index(lattice, 1, conj=false, forward=true, band=3) == 25
	@test index(lattice, 1, conj=true, forward=true, band=3) == 26
	@test index(lattice, 2, conj=false, forward=false, band=3) == 27
	@test index(lattice, 2, conj=true, forward=false, band=3) == 28
	@test index(lattice, 1, conj=false, forward=false, band=3) == 29
	@test index(lattice, 1, conj=true, forward=false, band=3) == 30
	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

@testset "GrassmannLattice: real time A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2" begin
	# one band
	lattice = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2())
	@test isa(ConjugationStyle(lattice), GeneralConjugation)
	@test isa(LayoutStyle(lattice), BranchLocalLayout)
	@test isa(lattice, RealGrassmannLattice)
	@test lattice.ordering == A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 14
	@test lattice.N == 2
	@test lattice.k == 3
	@test lattice.t == 0.1
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2

	@test index(lattice, 3, conj=false, forward=true) == 3
	@test index(lattice, 3, conj=true, forward=true) == 4
	@test index(lattice, 2, conj=false, forward=true) == 5
	@test index(lattice, 2, conj=true, forward=true) == 6
	@test index(lattice, 1, conj=false, forward=true) == 7
	@test index(lattice, 1, conj=true, forward=true) == 8

	@test index(lattice, 1, conj=false, forward=false) == 9
	@test index(lattice, 1, conj=true, forward=false) == 10
	@test index(lattice, 2, conj=false, forward=false) == 11
	@test index(lattice, 2, conj=true, forward=false) == 12
	@test index(lattice, 3, conj=false, forward=false) == 13
	@test index(lattice, 3, conj=true, forward=false) == 14

	# mps = GrassmannMPS(scalartype(lattice), length(lattice))
	# @test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# two bands
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=2, ordering=A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 20
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 2
	@test index(lattice, 0, conj=true, band=2) == 3
	@test index(lattice, 0, conj=true, band=1) == 4

	@test index(lattice, 2, conj=false, forward=true, band=1) == 5
	@test index(lattice, 2, conj=false, forward=true, band=2) == 6
	@test index(lattice, 2, conj=true, forward=true, band=2) == 7
	@test index(lattice, 2, conj=true, forward=true, band=1) == 8
	@test index(lattice, 1, conj=false, forward=true, band=1) == 9
	@test index(lattice, 1, conj=false, forward=true, band=2) == 10
	@test index(lattice, 1, conj=true, forward=true, band=2) == 11
	@test index(lattice, 1, conj=true, forward=true, band=1) == 12


	@test index(lattice, 1, conj=false, forward=false, band=1) == 13
	@test index(lattice, 1, conj=false, forward=false, band=2) == 14
	@test index(lattice, 1, conj=true, forward=false, band=2) == 15
	@test index(lattice, 1, conj=true, forward=false, band=1) == 16
	@test index(lattice, 2, conj=false, forward=false, band=1) == 17
	@test index(lattice, 2, conj=false, forward=false, band=2) == 18
	@test index(lattice, 2, conj=true, forward=false, band=2) == 19
	@test index(lattice, 2, conj=true, forward=false, band=1) == 20
	# mps = vacuumstate(lattice)
	# @test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# three bands
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=3, ordering=A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 30
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 2
	@test index(lattice, 0, conj=false, band=3) == 3
	@test index(lattice, 0, conj=true, band=3) == 4
	@test index(lattice, 0, conj=true, band=2) == 5
	@test index(lattice, 0, conj=true, band=1) == 6

	@test index(lattice, 2, conj=false, forward=true, band=1) == 7
	@test index(lattice, 2, conj=false, forward=true, band=2) == 8
	@test index(lattice, 2, conj=false, forward=true, band=3) == 9
	@test index(lattice, 2, conj=true, forward=true, band=3) == 10
	@test index(lattice, 2, conj=true, forward=true, band=2) == 11
	@test index(lattice, 2, conj=true, forward=true, band=1) == 12
	@test index(lattice, 1, conj=false, forward=true, band=1) == 13
	@test index(lattice, 1, conj=false, forward=true, band=2) == 14
	@test index(lattice, 1, conj=false, forward=true, band=3) == 15
	@test index(lattice, 1, conj=true, forward=true, band=3) == 16
	@test index(lattice, 1, conj=true, forward=true, band=2) == 17
	@test index(lattice, 1, conj=true, forward=true, band=1) == 18


	@test index(lattice, 1, conj=false, forward=false, band=1) == 19
	@test index(lattice, 1, conj=false, forward=false, band=2) == 20
	@test index(lattice, 1, conj=false, forward=false, band=3) == 21
	@test index(lattice, 1, conj=true, forward=false, band=3) == 22
	@test index(lattice, 1, conj=true, forward=false, band=2) == 23
	@test index(lattice, 1, conj=true, forward=false, band=1) == 24
	@test index(lattice, 2, conj=false, forward=false, band=1) == 25
	@test index(lattice, 2, conj=false, forward=false, band=2) == 26
	@test index(lattice, 2, conj=false, forward=false, band=3) == 27
	@test index(lattice, 2, conj=true, forward=false, band=3) == 28
	@test index(lattice, 2, conj=true, forward=false, band=2) == 29
	@test index(lattice, 2, conj=true, forward=false, band=1) == 30
	# mps = vacuumstate(lattice)
	# @test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end


@testset "GrassmannLattice: real time A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2" begin
	# one band
	lattice = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2())
	@test isa(ConjugationStyle(lattice), AdjacentConjugation)
	@test isa(LayoutStyle(lattice), BranchLocalLayout)
	@test isa(lattice, RealGrassmannLattice)
	@test lattice.ordering == A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 14
	@test lattice.N == 2
	@test lattice.k == 3
	@test lattice.t == 0.1
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2

	@test index(lattice, 3, conj=false, forward=true) == 3
	@test index(lattice, 3, conj=true, forward=true) == 4
	@test index(lattice, 2, conj=false, forward=true) == 5
	@test index(lattice, 2, conj=true, forward=true) == 6
	@test index(lattice, 1, conj=false, forward=true) == 7
	@test index(lattice, 1, conj=true, forward=true) == 8

	@test index(lattice, 1, conj=false, forward=false) == 9
	@test index(lattice, 1, conj=true, forward=false) == 10
	@test index(lattice, 2, conj=false, forward=false) == 11
	@test index(lattice, 2, conj=true, forward=false) == 12
	@test index(lattice, 3, conj=false, forward=false) == 13
	@test index(lattice, 3, conj=true, forward=false) == 14

	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# two bands
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=2, ordering=A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 20
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4

	@test index(lattice, 2, conj=false, forward=true, band=1) == 5
	@test index(lattice, 2, conj=true, forward=true, band=1) == 6
	@test index(lattice, 2, conj=false, forward=true, band=2) == 7
	@test index(lattice, 2, conj=true, forward=true, band=2) == 8
	@test index(lattice, 1, conj=false, forward=true, band=1) == 9
	@test index(lattice, 1, conj=true, forward=true, band=1) == 10
	@test index(lattice, 1, conj=false, forward=true, band=2) == 11
	@test index(lattice, 1, conj=true, forward=true, band=2) == 12

	@test index(lattice, 1, conj=false, forward=false, band=1) == 13
	@test index(lattice, 1, conj=true, forward=false, band=1) == 14
	@test index(lattice, 1, conj=false, forward=false, band=2) == 15
	@test index(lattice, 1, conj=true, forward=false, band=2) == 16
	@test index(lattice, 2, conj=false, forward=false, band=1) == 17
	@test index(lattice, 2, conj=true, forward=false, band=1) == 18
	@test index(lattice, 2, conj=false, forward=false, band=2) == 19
	@test index(lattice, 2, conj=true, forward=false, band=2) == 20

	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# three bands
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=3, ordering=A2A2B2B2A1A1B1B1a1a1b1b1a2a2b2b2())
	@test isa(lattice, RealGrassmannLattice)
	@test length(lattice) == 30
	@test lattice.N == 1
	@test lattice.k == 2
	@test lattice.t == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4
	@test index(lattice, 0, conj=false, band=3) == 5
	@test index(lattice, 0, conj=true, band=3) == 6

	@test index(lattice, 2, conj=false, forward=true, band=1) == 7
	@test index(lattice, 2, conj=true, forward=true, band=1) == 8
	@test index(lattice, 2, conj=false, forward=true, band=2) == 9
	@test index(lattice, 2, conj=true, forward=true, band=2) == 10
	@test index(lattice, 2, conj=false, forward=true, band=3) == 11
	@test index(lattice, 2, conj=true, forward=true, band=3) == 12
	@test index(lattice, 1, conj=false, forward=true, band=1) == 13
	@test index(lattice, 1, conj=true, forward=true, band=1) == 14
	@test index(lattice, 1, conj=false, forward=true, band=2) == 15
	@test index(lattice, 1, conj=true, forward=true, band=2) == 16
	@test index(lattice, 1, conj=false, forward=true, band=3) == 17
	@test index(lattice, 1, conj=true, forward=true, band=3) == 18

	@test index(lattice, 1, conj=false, forward=false, band=1) == 19
	@test index(lattice, 1, conj=true, forward=false, band=1) == 20
	@test index(lattice, 1, conj=false, forward=false, band=2) == 21
	@test index(lattice, 1, conj=true, forward=false, band=2) == 22
	@test index(lattice, 1, conj=false, forward=false, band=3) == 23
	@test index(lattice, 1, conj=true, forward=false, band=3) == 24
	@test index(lattice, 2, conj=false, forward=false, band=1) == 25
	@test index(lattice, 2, conj=true, forward=false, band=1) == 26
	@test index(lattice, 2, conj=false, forward=false, band=2) == 27
	@test index(lattice, 2, conj=true, forward=false, band=2) == 28
	@test index(lattice, 2, conj=false, forward=false, band=3) == 29
	@test index(lattice, 2, conj=true, forward=false, band=3) == 30
	
	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end
