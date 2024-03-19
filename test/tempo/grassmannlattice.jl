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
	@test index(lattice, 3, conj=false, branch=:+) == 3
	@test index(lattice, 3, conj=false, branch=:-) == 4
	@test index(lattice, 3, conj=true, branch=:+) == 6
	@test index(lattice, 3, conj=true, branch=:-) == 5
	@test index(lattice, 2, conj=false, branch=:+) == 7
	@test index(lattice, 2, conj=false, branch=:-) == 8
	@test index(lattice, 2, conj=true, branch=:+) == 10
	@test index(lattice, 2, conj=true, branch=:-) == 9	
	@test index(lattice, 1, conj=false, branch=:+) == 11
	@test index(lattice, 1, conj=false, branch=:-) == 12
	@test index(lattice, 1, conj=true, branch=:+) == 14
	@test index(lattice, 1, conj=true, branch=:-) == 13	

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 5
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 6
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 7
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 8
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 12
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 11
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 10
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 9

	@test index(lattice, 1, conj=false, branch=:+, band=1) == 13
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 15
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 16
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 20
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 19
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 18
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 17

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 7
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 9
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 10
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 11
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 12
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 18
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 17
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 16
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 15
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 14
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 13

	@test index(lattice, 1, conj=false, branch=:+, band=1) == 19
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 21
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 22
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 23
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 24
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 30
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 29
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 28
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 27
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 26
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 25

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

@testset "GrassmannLattice: real time A1B1ā1b̄1A1B1a1b1" begin
	# aā a₃^+ā₃^-ā₃^+a₃^- a₂^+ā₂^-ā₂^+a₂^-  a₁^+ā₁^-ā₁^+a₁^-
	lattice = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A1B1ā1b̄1A1B1a1b1())
	@test isa(ConjugationStyle(lattice), GeneralConjugation)
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, RealGrassmannLattice)
	@test lattice.ordering == A1B1ā1b̄1A1B1a1b1()
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
	@test index(lattice, 3, conj=false, branch=:+) == 3
	@test index(lattice, 3, conj=true, branch=:-) == 4
	@test index(lattice, 3, conj=true, branch=:+) == 5
	@test index(lattice, 3, conj=false, branch=:-) == 6
	
	@test index(lattice, 2, conj=false, branch=:+) == 7
	@test index(lattice, 2, conj=true, branch=:-) == 8
	@test index(lattice, 2, conj=true, branch=:+) == 9
	@test index(lattice, 2, conj=false, branch=:-) == 10
	@test index(lattice, 1, conj=false, branch=:+) == 11
	@test index(lattice, 1, conj=true, branch=:-) == 12
	@test index(lattice, 1, conj=true, branch=:+) == 13
	@test index(lattice, 1, conj=false, branch=:-) == 14

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# aābb̄ a₂^+b₂^+ā₂^-b̄₂^-ā₂^+b̄₂^+a₂^-b₂^-  a₁^+b₁^+ā₁^-b̄₁^-ā₁^+b̄₁^+a₁^-b₁^-
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=2, ordering=A1B1ā1b̄1A1B1a1b1())
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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 5
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 6
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 7
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 8
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 9
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 10
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 11
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 12

	@test index(lattice, 1, conj=false, branch=:+, band=1) == 13
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 14
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 15
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 16
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 17
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 18
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 19
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 20

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# aābb̄cc̄ a₂^+b₂^+c₂^+ā₂^-b̄₂^-c̄₂^-ā₂^+b̄₂^+c̄₂^+a₂^-b₂^-c₂^-  a₁^+b₁^+c₁^+ā₁^-b̄₁^-c̄₁^-ā₁^+b̄₁^+c̄₁^+a₁^-b₁^-c₁^-
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=3, ordering=A1B1ā1b̄1A1B1a1b1())
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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 7
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 8
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 9
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 10
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 12
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 13
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 14
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 15
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 16
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 17
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 18

	@test index(lattice, 1, conj=false, branch=:+, band=1) == 19
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 20
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 21
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 22
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 23
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 24
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 25
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 26
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 27
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 28
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 29
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 30

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
	@test index(lattice, 3, conj=false, branch=:+) == 3
	@test index(lattice, 3, conj=true, branch=:+) == 4
	@test index(lattice, 3, conj=false, branch=:-) == 5
	@test index(lattice, 3, conj=true, branch=:-) == 6
	@test index(lattice, 2, conj=false, branch=:+) == 7
	@test index(lattice, 2, conj=true, branch=:+) == 8
	@test index(lattice, 2, conj=false, branch=:-) == 9
	@test index(lattice, 2, conj=true, branch=:-) == 10
	@test index(lattice, 1, conj=false, branch=:+) == 11
	@test index(lattice, 1, conj=true, branch=:+) == 12
	@test index(lattice, 1, conj=false, branch=:-) == 13
	@test index(lattice, 1, conj=true, branch=:-) == 14

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 5
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 6
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 7
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 9
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 10
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 12

	@test index(lattice, 1, conj=false, branch=:+, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 15
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 16
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 17
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 18
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 19
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 20

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 7
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 9
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 10
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 12
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 13
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 14
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 15
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 16
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 17
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 18

	@test index(lattice, 1, conj=false, branch=:+, band=1) == 19
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 21
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 22
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 23
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 24
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 25
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 26
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 27
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 28
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 29
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 30

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

@testset "GrassmannLattice: real time A1A1B1B1a1a1b1b1" begin
	# aā a₁^+ā₁^+a₁^-ā₁^-
	lattice = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A1A1B1B1a1a1b1b1())
	@test isa(ConjugationStyle(lattice), AdjacentConjugation)
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, RealGrassmannLattice)
	@test lattice.ordering == A1A1B1B1a1a1b1b1()
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
	@test index(lattice, 3, conj=false, branch=:+) == 3
	@test index(lattice, 3, conj=true, branch=:+) == 4
	@test index(lattice, 3, conj=false, branch=:-) == 5
	@test index(lattice, 3, conj=true, branch=:-) == 6
	@test index(lattice, 2, conj=false, branch=:+) == 7
	@test index(lattice, 2, conj=true, branch=:+) == 8
	@test index(lattice, 2, conj=false, branch=:-) == 9
	@test index(lattice, 2, conj=true, branch=:-) == 10
	@test index(lattice, 1, conj=false, branch=:+) == 11
	@test index(lattice, 1, conj=true, branch=:+) == 12
	@test index(lattice, 1, conj=false, branch=:-) == 13
	@test index(lattice, 1, conj=true, branch=:-) == 14

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# aābb̄ a₂^+ā₂^+b₂^+b̄₂^+a₂^-ā₂^-b₂^-b̄₂^- a₁^+ā₁^+b₁^+b̄₁^+a₁^-ā₁^-b₁^-b̄₁^-
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=2, ordering=A1A1B1B1a1a1b1b1())
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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 5
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 6
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 7
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 8
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 9
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 10
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 12

	@test index(lattice, 1, conj=false, branch=:+, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 16
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 17
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 18
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 19
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 20

	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# aābb̄cc̄ a₂^+ā₂^+b₂^+b̄₂^+c₂^+c̄₂^+a₂^-ā₂^-b₂^-b̄₂^-c₂^-c̄₂^- a₁^+ā₁^+b₁^+b̄₁^+c₁^+c̄₁^+a₁^-ā₁^-b₁^-b̄₁^-c₁^-c̄₁^-
	lattice = GrassmannLattice(N=1, δt=0.1, contour=:real, bands=3, ordering=A1A1B1B1a1a1b1b1())
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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 7
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 9
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 10
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 11
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 12
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 13
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 14
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 15
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 16
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 17
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 18

	@test index(lattice, 1, conj=false, branch=:+, band=1) == 19
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 21
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 22
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 23
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 24
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 25
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 26
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 27
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 28
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 29
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 30

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
	@test index(lattice, 3, conj=false, branch=:+) == 3
	@test index(lattice, 3, conj=true, branch=:+) == 4
	@test index(lattice, 2, conj=false, branch=:+) == 5
	@test index(lattice, 2, conj=true, branch=:+) == 6
	@test index(lattice, 1, conj=false, branch=:+) == 7
	@test index(lattice, 1, conj=true, branch=:+) == 8
	@test index(lattice, 3, conj=false, branch=:-) == 9
	@test index(lattice, 3, conj=true, branch=:-) == 10
	@test index(lattice, 2, conj=false, branch=:-) == 11
	@test index(lattice, 2, conj=true, branch=:-) == 12
	@test index(lattice, 1, conj=false, branch=:-) == 13
	@test index(lattice, 1, conj=true, branch=:-) == 14

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 5
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 6
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 7
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 9
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 10
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 11
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 12

	@test index(lattice, 2, conj=false, branch=:+, band=2) == 13
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 14
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 16
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 17
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 18
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 19
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 20

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 7
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 8
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 9
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 10
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 11
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 12
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 14

	@test index(lattice, 2, conj=false, branch=:+, band=2) == 15
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 16
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 17
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 18
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 19
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 20
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 21
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 22

	@test index(lattice, 2, conj=false, branch=:+, band=3) == 23
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 24
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 25
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 26
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 27
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 28
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 29
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 30
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

	@test index(lattice, 3, conj=false, branch=:+) == 3
	@test index(lattice, 3, conj=true, branch=:+) == 4
	@test index(lattice, 2, conj=false, branch=:+) == 5
	@test index(lattice, 2, conj=true, branch=:+) == 6
	@test index(lattice, 1, conj=false, branch=:+) == 7
	@test index(lattice, 1, conj=true, branch=:+) == 8

	@test index(lattice, 1, conj=false, branch=:-) == 9
	@test index(lattice, 1, conj=true, branch=:-) == 10
	@test index(lattice, 2, conj=false, branch=:-) == 11
	@test index(lattice, 2, conj=true, branch=:-) == 12
	@test index(lattice, 3, conj=false, branch=:-) == 13
	@test index(lattice, 3, conj=true, branch=:-) == 14

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 5
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 6
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 7
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 8
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 9
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 10
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 11
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 12


	@test index(lattice, 1, conj=false, branch=:-, band=1) == 13
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 14
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 16
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 17
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 18
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 19
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 20
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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 7
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 8
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 9
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 10
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 12
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 13
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 14
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 15
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 16
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 17
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 18


	@test index(lattice, 1, conj=false, branch=:-, band=1) == 19
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 20
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 21
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 22
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 23
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 24
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 25
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 26
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 27
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 28
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 29
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 30
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

	@test index(lattice, 3, conj=false, branch=:+) == 3
	@test index(lattice, 3, conj=true, branch=:+) == 4
	@test index(lattice, 2, conj=false, branch=:+) == 5
	@test index(lattice, 2, conj=true, branch=:+) == 6
	@test index(lattice, 1, conj=false, branch=:+) == 7
	@test index(lattice, 1, conj=true, branch=:+) == 8

	@test index(lattice, 1, conj=false, branch=:-) == 9
	@test index(lattice, 1, conj=true, branch=:-) == 10
	@test index(lattice, 2, conj=false, branch=:-) == 11
	@test index(lattice, 2, conj=true, branch=:-) == 12
	@test index(lattice, 3, conj=false, branch=:-) == 13
	@test index(lattice, 3, conj=true, branch=:-) == 14

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 5
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 6
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 7
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 8
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 9
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 10
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 11
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 12

	@test index(lattice, 1, conj=false, branch=:-, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 16
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 17
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 18
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 19
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 20

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

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 7
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 9
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 10
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 11
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 12
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 16
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 17
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 18

	@test index(lattice, 1, conj=false, branch=:-, band=1) == 19
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 21
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 22
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 23
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 24
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 25
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 26
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 27
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 28
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 29
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 30
	
	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

@testset "GrassmannLattice: mixed time A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2" begin
	# imag part: a\bar{a} a_2\bar{a}_2 a_1\bar{a}_1
	# real part: a_3^+\bar{a}_3^+a_3^-\bar{a}_3^- a_2^+\bar{a}_2^+a_2^-\bar{a}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-
	lattice = GrassmannLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2())
	@test isa(ConjugationStyle(lattice), AdjacentConjugation)
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, MixedGrassmannLattice)
	@test lattice.ordering == A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 20
	@test lattice.Nt == 2
	@test lattice.Nτ == 2
	@test lattice.t == 0.1
	@test lattice.β == 0.2
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test lattice.τs == 0:0.1:0.2

	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2

	@test index(lattice, 3, conj=false, branch=:τ) == 3
	@test index(lattice, 3, conj=true, branch=:τ) == 4
	@test index(lattice, 2, conj=false, branch=:τ) == 5
	@test index(lattice, 2, conj=true, branch=:τ) == 6	
	@test index(lattice, 1, conj=false, branch=:τ) == 7 
	@test index(lattice, 1, conj=true, branch=:τ) == 8


	@test index(lattice, 1, conj=false, branch=:+) == 9
	@test index(lattice, 1, conj=true, branch=:+) == 10
	@test index(lattice, 1, conj=false, branch=:-) == 11
	@test index(lattice, 1, conj=true, branch=:-) == 12
	@test index(lattice, 2, conj=false, branch=:+) == 13
	@test index(lattice, 2, conj=true, branch=:+) == 14
	@test index(lattice, 2, conj=false, branch=:-) == 15
	@test index(lattice, 2, conj=true, branch=:-) == 16
	@test index(lattice, 3, conj=false, branch=:+) == 17
	@test index(lattice, 3, conj=true, branch=:+) == 18
	@test index(lattice, 3, conj=false, branch=:-) == 19
	@test index(lattice, 3, conj=true, branch=:-) == 20


	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# imag part: a\bar{a}b\bar{b} a_2\bar{a}_2b_2\bar{b}_2 a_1\bar{a}_1b_1\bar{b}_1
	# real part: a_2^+\bar{a}_2^+a_2^-\bar{a}_2^-b_2^+\bar{b}_2^+b_2^-\bar{b}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-b_1^+\bar{b}_1^+b_1^-\bar{b}_1^-
	lattice = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, contour=:mixed, bands=2, ordering=A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2())
	@test isa(lattice, MixedGrassmannLattice)
	@test length(lattice) == 32
	@test lattice.Nt == 1
	@test lattice.Nτ == 2
	@test lattice.t == 0.1
	@test lattice.β == 0.2
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test lattice.τs == 0:0.1:0.2

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4

	@test index(lattice, 3, conj=false, branch=:τ, band=1) == 5
	@test index(lattice, 3, conj=true, branch=:τ, band=1) == 6
	@test index(lattice, 3, conj=false, branch=:τ, band=2) == 7
	@test index(lattice, 3, conj=true, branch=:τ, band=2) == 8
	@test index(lattice, 2, conj=false, branch=:τ, band=1) == 9
	@test index(lattice, 2, conj=true, branch=:τ, band=1) == 10
	@test index(lattice, 2, conj=false, branch=:τ, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:τ, band=2) == 12
	@test index(lattice, 1, conj=false, branch=:τ, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:τ, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:τ, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:τ, band=2) == 16


	@test index(lattice, 1, conj=false, branch=:+, band=1) == 17
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 18
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 19
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 21
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 22
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 23
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 24

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 25
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 26
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 27
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 28
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 29
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 30
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 31
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 32

	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# imag part: a\bar{a}b\bar{b}c\bar{c} a_1\bar{a}_1b_1\bar{b}_1c_1\bar{c}_1
	# real part: a_2^+\bar{a}_2^+a_2^-\bar{a}_2^-b_2^+\bar{b}_2^+b_2^-\bar{b}_2^-c_2^+\bar{c}_2^+c_2^-\bar{c}_2^- 
	# a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-b_1^+\bar{b}_1^+b_1^-\bar{b}_1^-c_1^+\bar{c}_1^+c_1^-\bar{c}_1^-
	lattice = GrassmannLattice(Nt=1, δt=0.1, Nτ=1, δτ=0.1, contour=:mixed, bands=3, ordering=A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2())
	@test isa(lattice, MixedGrassmannLattice)
	@test length(lattice) == 42
	@test lattice.Nt == 1
	@test lattice.Nτ == 1
	@test lattice.t == 0.1
	@test lattice.β == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test lattice.τs == 0:0.1:0.1

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4
	@test index(lattice, 0, conj=false, band=3) == 5
	@test index(lattice, 0, conj=true, band=3) == 6

	@test index(lattice, 2, conj=false, branch=:τ, band=1) == 7
	@test index(lattice, 2, conj=true, branch=:τ, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:τ, band=2) == 9
	@test index(lattice, 2, conj=true, branch=:τ, band=2) == 10
	@test index(lattice, 2, conj=false, branch=:τ, band=3) == 11
	@test index(lattice, 2, conj=true, branch=:τ, band=3) == 12
	@test index(lattice, 1, conj=false, branch=:τ, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:τ, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:τ, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:τ, band=2) == 16
	@test index(lattice, 1, conj=false, branch=:τ, band=3) == 17
	@test index(lattice, 1, conj=true, branch=:τ, band=3) == 18


	@test index(lattice, 1, conj=false, branch=:+, band=1) == 19
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:-, band=1) == 21
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 22
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 23
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 24
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 25
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 26
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 27
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 28
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 29
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 30

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 31
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 32
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 33
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 34
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 35
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 36
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 37
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 38
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 39
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 40
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 41
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 42

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

# swap forward and backward branches compared to A1A1B1B1_A1A1a1a1B1B1b1b1A2A2a2a2B2B2b2b2 
@testset "GrassmannLattice: mixed time A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2" begin
	# one band
	lattice = GrassmannLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2())
	@test isa(ConjugationStyle(lattice), AdjacentConjugation)
	@test isa(LayoutStyle(lattice), TimeLocalLayout)
	@test isa(lattice, MixedGrassmannLattice)
	@test lattice.ordering == A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 20
	@test lattice.Nt == 2
	@test lattice.Nτ == 2
	@test lattice.t == 0.1
	@test lattice.β == 0.2
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test lattice.τs == 0:0.1:0.2

	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2

	@test index(lattice, 3, conj=false, branch=:τ) == 3
	@test index(lattice, 3, conj=true, branch=:τ) == 4
	@test index(lattice, 2, conj=false, branch=:τ) == 5
	@test index(lattice, 2, conj=true, branch=:τ) == 6	
	@test index(lattice, 1, conj=false, branch=:τ) == 7
	@test index(lattice, 1, conj=true, branch=:τ) == 8

	@test index(lattice, 1, conj=false, branch=:-) == 9
	@test index(lattice, 1, conj=true, branch=:-) == 10
	@test index(lattice, 1, conj=false, branch=:+) == 11
	@test index(lattice, 1, conj=true, branch=:+) == 12
	@test index(lattice, 2, conj=false, branch=:-) == 13
	@test index(lattice, 2, conj=true, branch=:-) == 14
	@test index(lattice, 2, conj=false, branch=:+) == 15
	@test index(lattice, 2, conj=true, branch=:+) == 16
	@test index(lattice, 3, conj=false, branch=:-) == 17
	@test index(lattice, 3, conj=true, branch=:-) == 18
	@test index(lattice, 3, conj=false, branch=:+) == 19
	@test index(lattice, 3, conj=true, branch=:+) == 20

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# two bands
	lattice = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, contour=:mixed, bands=2, ordering=A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2())
	@test isa(lattice, MixedGrassmannLattice)
	@test length(lattice) == 32
	@test lattice.Nt == 1
	@test lattice.Nτ == 2
	@test lattice.t == 0.1
	@test lattice.β == 0.2
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test lattice.τs == 0:0.1:0.2

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4

	@test index(lattice, 3, conj=false, branch=:τ, band=1) == 5
	@test index(lattice, 3, conj=true, branch=:τ, band=1) == 6
	@test index(lattice, 3, conj=false, branch=:τ, band=2) == 7
	@test index(lattice, 3, conj=true, branch=:τ, band=2) == 8
	@test index(lattice, 2, conj=false, branch=:τ, band=1) == 9
	@test index(lattice, 2, conj=true, branch=:τ, band=1) == 10
	@test index(lattice, 2, conj=false, branch=:τ, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:τ, band=2) == 12
	@test index(lattice, 1, conj=false, branch=:τ, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:τ, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:τ, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:τ, band=2) == 16

	@test index(lattice, 1, conj=false, branch=:-, band=1) == 17
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 18
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 19
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 21
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 22
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 23
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 24

	@test index(lattice, 2, conj=false, branch=:-, band=1) == 25
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 26
	@test index(lattice, 2, conj=false, branch=:+, band=1) == 27
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 28
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 29
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 30
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 31
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 32

	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# three bands
	lattice = GrassmannLattice(Nt=1, δt=0.1, Nτ=1, δτ=0.1, contour=:mixed, bands=3, ordering=A1A1B1B1_a1a1A1A1b1b1B1B1a2a2A2A2b2b2B2B2())
	@test isa(lattice, MixedGrassmannLattice)
	@test length(lattice) == 42
	@test lattice.Nt == 1
	@test lattice.Nτ == 1
	@test lattice.t == 0.1
	@test lattice.β == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test lattice.τs == 0:0.1:0.1

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=true, band=1) == 2
	@test index(lattice, 0, conj=false, band=2) == 3
	@test index(lattice, 0, conj=true, band=2) == 4
	@test index(lattice, 0, conj=false, band=3) == 5
	@test index(lattice, 0, conj=true, band=3) == 6

	@test index(lattice, 2, conj=false, branch=:τ, band=1) == 7
	@test index(lattice, 2, conj=true, branch=:τ, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:τ, band=2) == 9
	@test index(lattice, 2, conj=true, branch=:τ, band=2) == 10
	@test index(lattice, 2, conj=false, branch=:τ, band=3) == 11
	@test index(lattice, 2, conj=true, branch=:τ, band=3) == 12
	@test index(lattice, 1, conj=false, branch=:τ, band=1) == 13
	@test index(lattice, 1, conj=true, branch=:τ, band=1) == 14
	@test index(lattice, 1, conj=false, branch=:τ, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:τ, band=2) == 16
	@test index(lattice, 1, conj=false, branch=:τ, band=3) == 17
	@test index(lattice, 1, conj=true, branch=:τ, band=3) == 18

	@test index(lattice, 1, conj=false, branch=:-, band=1) == 19
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 21
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 22
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 23
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 24
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 25
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 26
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 27
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 28
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 29
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 30

	@test index(lattice, 2, conj=false, branch=:-, band=1) == 31
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 32
	@test index(lattice, 2, conj=false, branch=:+, band=1) == 33
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 34
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 35
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 36
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 37
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 38
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 39
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 40
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 41
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 42

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

# mixed time lattice
@testset "GrassmannLattice: mixed time A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2" begin
	# imag part: a_2\bar{a}_2 a_1\bar{a}_1
	# real part: a\bar{a} a_3^+\bar{a}_3^+a_3^-\bar{a}_3^- a_2^+\bar{a}_2^+a_2^-\bar{a}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-
	lattice = GrassmannLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2())
	@test isa(ConjugationStyle(lattice), GeneralConjugation)
	@test isa(LayoutStyle(lattice), BranchLocalLayout)
	@test isa(lattice, MixedGrassmannLattice)
	@test lattice.ordering == A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2()
	@test scalartype(lattice) == ComplexF64
	@test length(lattice) == 20
	@test lattice.Nt == 2
	@test lattice.Nτ == 2
	@test lattice.t == 0.1
	@test lattice.β == 0.2
	@test lattice.bands == 1
	@test lattice.δt == 0.05
	@test lattice.ts == 0:0.05:0.1
	@test lattice.τs == 0:0.1:0.2

	@test index(lattice, 0, conj=false) == 1
	@test index(lattice, 0, conj=true) == 2

	@test index(lattice, 3, conj=false, branch=:τ) == 3
	@test index(lattice, 3, conj=true, branch=:τ) == 4
	@test index(lattice, 2, conj=false, branch=:τ) == 5
	@test index(lattice, 2, conj=true, branch=:τ) == 6
	@test index(lattice, 1, conj=false, branch=:τ) == 7
	@test index(lattice, 1, conj=true, branch=:τ) == 8

	@test index(lattice, 3, conj=false, branch=:+) == 9
	@test index(lattice, 3, conj=true, branch=:+) == 10
	@test index(lattice, 2, conj=false, branch=:+) == 11
	@test index(lattice, 2, conj=true, branch=:+) == 12
	@test index(lattice, 1, conj=false, branch=:+) == 13
	@test index(lattice, 1, conj=true, branch=:+) == 14

	@test index(lattice, 1, conj=false, branch=:-) == 15
	@test index(lattice, 1, conj=true, branch=:-) == 16
	@test index(lattice, 2, conj=false, branch=:-) == 17
	@test index(lattice, 2, conj=true, branch=:-) == 18
	@test index(lattice, 3, conj=false, branch=:-) == 19
	@test index(lattice, 3, conj=true, branch=:-) == 20

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# imag part: a_2b_2\bar{b}_2\bar{a}_2 a_1b_1\bar{b}_1\bar{a}_1
	# real part: a\bar{a}b\bar{b} a_2^+\bar{a}_2^+a_2^-\bar{a}_2^-b_2^+\bar{b}_2^+b_2^-\bar{b}_2^- a_1^+\bar{a}_1^+a_1^-\bar{a}_1^-b_1^+\bar{b}_1^+b_1^-\bar{b}_1^-
	lattice = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, contour=:mixed, bands=2, ordering=A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2())
	@test isa(lattice, MixedGrassmannLattice)
	@test length(lattice) == 32
	@test lattice.Nt == 1
	@test lattice.Nτ == 2
	@test lattice.t == 0.1
	@test lattice.β == 0.2
	@test lattice.bands == 2
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test lattice.τs == 0:0.1:0.2

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 2
	@test index(lattice, 0, conj=true, band=2) == 3
	@test index(lattice, 0, conj=true, band=1) == 4

	@test index(lattice, 3, conj=false, branch=:τ, band=1) == 5
	@test index(lattice, 3, conj=false, branch=:τ, band=2) == 6
	@test index(lattice, 3, conj=true, branch=:τ, band=2) == 7
	@test index(lattice, 3, conj=true, branch=:τ, band=1) == 8
	@test index(lattice, 2, conj=false, branch=:τ, band=1) == 9
	@test index(lattice, 2, conj=false, branch=:τ, band=2) == 10
	@test index(lattice, 2, conj=true, branch=:τ, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:τ, band=1) == 12
	@test index(lattice, 1, conj=false, branch=:τ, band=1) == 13
	@test index(lattice, 1, conj=false, branch=:τ, band=2) == 14
	@test index(lattice, 1, conj=true, branch=:τ, band=2) == 15
	@test index(lattice, 1, conj=true, branch=:τ, band=1) == 16

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 17
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 18
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 19
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 20
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 21
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 22
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 23
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 24

	@test index(lattice, 1, conj=false, branch=:-, band=1) == 25
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 26
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 27
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 28
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 29
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 30
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 31
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 32

	mps = vacuumstate(lattice)
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6

	# imag part: a_1b_1c_1\bar{c}_1\bar{b}_1\bar{a}_1
	lattice = GrassmannLattice(Nt=1, δt=0.1, Nτ=1, δτ=0.1, contour=:mixed, bands=3, ordering=A1B1B1A1_A2B2B2A2A1B1B1A1a1b1b1a1a2b2b2a2())
	@test isa(lattice, MixedGrassmannLattice)
	@test length(lattice) == 42
	@test lattice.Nt == 1
	@test lattice.Nτ == 1
	@test lattice.t == 0.1
	@test lattice.β == 0.1
	@test lattice.bands == 3
	@test lattice.δt == 0.1
	@test lattice.ts == 0:0.1:0.1
	@test lattice.τs == 0:0.1:0.1

	@test index(lattice, 0, conj=false, band=1) == 1
	@test index(lattice, 0, conj=false, band=2) == 2
	@test index(lattice, 0, conj=false, band=3) == 3
	@test index(lattice, 0, conj=true, band=3) == 4
	@test index(lattice, 0, conj=true, band=2) == 5
	@test index(lattice, 0, conj=true, band=1) == 6

	@test index(lattice, 2, conj=false, branch=:τ, band=1) == 7
	@test index(lattice, 2, conj=false, branch=:τ, band=2) == 8
	@test index(lattice, 2, conj=false, branch=:τ, band=3) == 9
	@test index(lattice, 2, conj=true, branch=:τ, band=3) == 10
	@test index(lattice, 2, conj=true, branch=:τ, band=2) == 11
	@test index(lattice, 2, conj=true, branch=:τ, band=1) == 12
	@test index(lattice, 1, conj=false, branch=:τ, band=1) == 13
	@test index(lattice, 1, conj=false, branch=:τ, band=2) == 14
	@test index(lattice, 1, conj=false, branch=:τ, band=3) == 15
	@test index(lattice, 1, conj=true, branch=:τ, band=3) == 16
	@test index(lattice, 1, conj=true, branch=:τ, band=2) == 17
	@test index(lattice, 1, conj=true, branch=:τ, band=1) == 18

	@test index(lattice, 2, conj=false, branch=:+, band=1) == 19
	@test index(lattice, 2, conj=false, branch=:+, band=2) == 20
	@test index(lattice, 2, conj=false, branch=:+, band=3) == 21
	@test index(lattice, 2, conj=true, branch=:+, band=3) == 22
	@test index(lattice, 2, conj=true, branch=:+, band=2) == 23
	@test index(lattice, 2, conj=true, branch=:+, band=1) == 24
	@test index(lattice, 1, conj=false, branch=:+, band=1) == 25
	@test index(lattice, 1, conj=false, branch=:+, band=2) == 26
	@test index(lattice, 1, conj=false, branch=:+, band=3) == 27
	@test index(lattice, 1, conj=true, branch=:+, band=3) == 28
	@test index(lattice, 1, conj=true, branch=:+, band=2) == 29
	@test index(lattice, 1, conj=true, branch=:+, band=1) == 30

	@test index(lattice, 1, conj=false, branch=:-, band=1) == 31
	@test index(lattice, 1, conj=false, branch=:-, band=2) == 32
	@test index(lattice, 1, conj=false, branch=:-, band=3) == 33
	@test index(lattice, 1, conj=true, branch=:-, band=3) == 34
	@test index(lattice, 1, conj=true, branch=:-, band=2) == 35
	@test index(lattice, 1, conj=true, branch=:-, band=1) == 36
	@test index(lattice, 2, conj=false, branch=:-, band=1) == 37
	@test index(lattice, 2, conj=false, branch=:-, band=2) == 38
	@test index(lattice, 2, conj=false, branch=:-, band=3) == 39
	@test index(lattice, 2, conj=true, branch=:-, band=3) == 40
	@test index(lattice, 2, conj=true, branch=:-, band=2) == 41
	@test index(lattice, 2, conj=true, branch=:-, band=1) == 42

	mps = GrassmannMPS(scalartype(lattice), length(lattice))
	@test integrate(lattice, mps) ≈ 1 atol = 1.0e-6
end

@testset "ImagGrassmannLattice: change ordering" begin
	function check_lattice_match(tsc, src, mapping)
		for i in 0:tsc.k
			for c in (true, false)
				for band in 1:tsc.bands
					if mapping[index(tsc, i, conj=c, band=band)] != index(src, i, conj=c, band=band)
						return false
					end
				end
			end
		end
		return true
	end
	
	N = 3
	bands = 2
	for o1 in imag_grassmann_orderings
		tsc = GrassmannLattice(δτ=0.1, N=N, contour=:imag, bands=bands, ordering=o1)
		for o2 in imag_grassmann_orderings
			src = GrassmannLattice(δτ=0.1, N=N, contour=:imag, bands=bands, ordering=o2)
			mapping = matchindices(tsc, src)
			@test check_lattice_match(tsc, src, mapping)
		end
	end
end


@testset "RealGrassmannLattice: change ordering" begin
	function check_lattice_match(tsc, src, mapping)
		for i in 0:tsc.k
			for c in (true, false)
				for band in 1:tsc.bands
					for f in (:+, :-)
						if mapping[index(tsc, i, conj=c, branch=f, band=band)] != index(src, i, conj=c, branch=f, band=band)
							return false
						end
					end
				end
			end
		end
		return true
	end
	
	N = 2
	bands = 3
	for o1 in real_grassmann_orderings
		tsc = GrassmannLattice(δt=0.1, N=N, contour=:real, bands=bands, ordering=o1)
		for o2 in real_grassmann_orderings
			src = GrassmannLattice(δt=0.1, N=N, contour=:real, bands=bands, ordering=o2)
			mapping = matchindices(tsc, src)
			@test check_lattice_match(tsc, src, mapping)
		end
	end
end


@testset "MixedGrassmannLattice: change ordering" begin
	function check_lattice_match(tsc, src, mapping)
		for i in 1:tsc.Nτ
			for c in (true, false)
				for band in 1:tsc.bands
					f = :τ
					if mapping[index(tsc, i, conj=c, branch=f, band=band)] != index(src, i, conj=c, branch=f, band=band)
						return false
					end
				end
			end
		end
		for i in 0:tsc.Nt+1
			for c in (true, false)
				for band in 1:tsc.bands
					for f in (:+, :-)
						if mapping[index(tsc, i, conj=c, branch=f, band=band)] != index(src, i, conj=c, branch=f, band=band)
							return false
						end
					end
				end
			end
		end
		return true
	end
	
	Nt = 1
	Nτ = 2
	bands = 3
	for o1 in mixed_grassmann_orderings
		tsc = GrassmannLattice(δt=0.1, Nt=Nt, δτ=0.03, Nτ=Nτ, contour=:mixed, bands=bands, ordering=o1)
		for o2 in mixed_grassmann_orderings
			src = similar(tsc, ordering=o2)
			mapping = matchindices(tsc, src)
			@test check_lattice_match(tsc, src, mapping)
		end
	end
end

@testset "GrassmannLattice: swapbandperm" begin
	N = 2
	bands = 3

	change_band(b, b1, b2) = ifelse(b == b1, b2, ifelse(b == b2, b1, b))

	for o1 in imag_grassmann_orderings
		x = GrassmannLattice(δτ=0.1, N=N, contour=:imag, bands=bands, ordering=o1)
		perm = swapbandperm(x, 1, 2)
		m = indexmappings(x)
		for ((j, c, b, band), pos) in m
			@test m[(j, c, b, change_band(band, 1, 2))] == perm[pos]
		end

		# mps = randomgmps(scalartype(x), length(x), D=4)
		# mps2 = swapband(swapband(mps, x, 1, 2), x, 1, 2)
		# @test distance(mps2, mps) / norm(mps) < 1.0e-6
	end

	for o1 in real_grassmann_orderings
		x = GrassmannLattice(δt=0.1, N=N, contour=:real, bands=bands, ordering=o1)

		perm = swapbandperm(x, 1, 3)
		m = indexmappings(x)
		for ((j, c, b, band), pos) in m
			@test m[(j, c, b, change_band(band, 1, 3))] == perm[pos]
		end		

		# mps = randomgmps(scalartype(x), length(x), D=4)
		# mps2 = swapband(swapband(mps, x, 1, 3), x, 1, 3)
		# @test distance(mps2, mps) / norm(mps) < 1.0e-6
	end

	Nt = 1
	Nτ = 2
	bands = 3
	for o1 in mixed_grassmann_orderings
		x = GrassmannLattice(δt=0.1, Nt=Nt, δτ=0.03, Nτ=Nτ, contour=:mixed, bands=bands, ordering=o1)

		perm = swapbandperm(x, 2, 3)
		m = indexmappings(x)
		for ((j, c, b, band), pos) in m
			@test m[(j, c, b, change_band(band, 3, 2))] == perm[pos]
		end		

		# mps = randomgmps(scalartype(x), length(x), D=4)
		# mps2 = swapband(swapband(mps, x, 2, 3), x, 3, 2)
		# @test distance(mps2, mps) / norm(mps) < 1.0e-6
	end
end

