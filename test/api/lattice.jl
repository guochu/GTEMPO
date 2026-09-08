@testset "API: GrassmannLattice" begin
	# imaginary-time contour
	lat = GrassmannLattice(N=6, δτ=0.1, contour=:imag)
	@test lat isa ImagGrassmannLattice1Order
	@test lat.k == lat.kτ == 7
	@test lat.N == lat.Nτ == 6
	@test lat.β ≈ 0.6
	@test collect(lat.τs) ≈ collect(0:0.1:0.6)
	@test lat.bands == 1
	@test length(lat) == 2 * lat.bands * (lat.k + 1)
	@test scalartype(lat) == Float64
	@test branches(lat) == (:τ,)

	lat2 = similar(lat, bands=2)
	@test lat2.bands == 2 && lat2.δτ == lat.δτ && lat2.N == lat.N
	lat3 = similar(lat, N=8, δτ=0.05)
	@test lat3.N == 8 && lat3.δτ == 0.05

	# real-time contour
	lat = GrassmannLattice(N=5, δt=0.1, contour=:real)
	@test lat isa RealGrassmannLattice1Order
	@test lat.k == lat.kt == 6
	@test lat.t == 0.5 && lat.Nt == 5
	@test lat.ts == 0:0.1:0.5
	@test scalartype(lat) == ComplexF64
	@test branches(lat) == (:+, :-)
	@test length(lat) == 4 * lat.bands * lat.k + 2 * lat.bands

	lat = GrassmannLattice(N=5, δt=0.1, contour=:Keldysh)
	@test lat isa RealGrassmannLattice1Order
	lat = GrassmannLattice(N=5, δt=0.1, contour=:real, order=2)
	@test lat isa RealGrassmannLattice2Order

	# mixed (Kadanoff) contour
	lat = GrassmannLattice(Nt=4, δt=0.1, Nτ=5, δτ=0.05, contour=:mixed)
	@test lat.kt == 5 && lat.kτ == 6
	@test lat.t ≈ 0.4 && lat.β ≈ 0.25
	@test branches(lat) == (:+, :-, :τ)
	@test length(lat) == 4*lat.bands*lat.kt + 2*lat.bands + 2*lat.bands*(lat.Nτ+1)

	lat = GrassmannLattice(Nt=4, δt=0.1, Nτ=5, δτ=0.05, contour=:Kadanoff)
	@test lat.kt == 5 && lat.kτ == 6

	# orderings
	for o in imag_orderings
		@test GrassmannLattice(N=4, δτ=0.1, contour=:imag, ordering=o).ordering === o
	end
	for o in real_orderings
		@test GrassmannLattice(N=4, δt=0.1, contour=:real, ordering=o).ordering === o
	end
	for o in mixed_orderings
		@test GrassmannLattice(Nt=3, δt=0.1, Nτ=4, δτ=0.1, contour=:mixed, ordering=o).ordering === o
	end

	# orbitals (bcs convention: bands = 2 * orbitals)
	@test GrassmannLattice(N=4, δτ=0.1, bands=2, contour=:imag).orbitals == 1
	@test GrassmannLattice(Nt=3, δt=0.1, Nτ=3, δτ=0.1, bands=4, contour=:mixed).orbitals == 2

	# index consistency: indexmappings must be a bijection over all variables
	lat = GrassmannLattice(N=3, δτ=0.1, bands=2, contour=:imag)
	m = indexmappings(lat)
	@test length(m) == (lat.k + 1) * 2 * lat.bands
	@test length(unique(values(m))) == length(m)

	lat = GrassmannLattice(N=3, δt=0.1, contour=:real, bands=2)
	m = indexmappings(lat)
	@test length(m) == 2 * lat.bands + 4 * lat.bands * lat.k  # i=0 has no branch
	@test length(unique(values(m))) == length(m)

	# ContourIndex: construction and contour ordering semantics
	ci1 = ContourIndex(2, conj=false, branch=:+, band=1)
	ci2 = ContourIndex(1, conj=true, branch=:+, band=1)
	ci3 = ContourIndex(3, conj=false, branch=:-, band=1)
	ci4 = ContourIndex(1, conj=false, branch=:τ, band=1)
	@test ci1 == ContourIndex(2, band=1, conj=false, branch=:+)
	@test ci2 < ci1                      # + branch: ascending time, conj first
	@test ci1 < ci3                      # + precedes -
	@test ci3 < ci4                      # - precedes τ
	@test ci4 < ContourIndex(3, conj=true, branch=:τ, band=1)
	@test branch(ci1) == :+
	@test_throws ArgumentError ci1 < ContourIndex(0, conj=true, branch=:+, band=1)
end

# all exported orderings of each contour, including the retarded-interaction
# orderings and the fast-propagator mixed-time ordering
const all_imag_orderings = [A1Ā1B1B̄1(), A1B1B̄1Ā1(), A2Ā2A1Ā1B2B̄2B1B̄1(), Ā2A1B̄2B1()]
const all_real_orderings = [A1Ā1B1B̄1a1ā1b1b̄1(), A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(),
							A1B1ā1b̄1Ā1B̄1a1b1(), A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1(),
							A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2(),
							Ā2A1ā1a2B̄2B1b̄1b̄2()]
const all_mixed_orderings = [A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2(),
							 A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2(),
							 A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(),
							 Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2(),
							 A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1()]

# ConjugationStyle/LayoutStyle classification of every exported ordering
# (ground truth: src/lattices/grassmannordering.jl)
@testset "API: ordering classification" begin
	for (o, conj, layout) in ((A1Ā1B1B̄1(), AdjacentConjugation, TimeLocalLayout),
							 (A1B1B̄1Ā1(), GeneralConjugation, TimeLocalLayout),
							 (A2Ā2A1Ā1B2B̄2B1B̄1(), AdjacentConjugation, BandLocalLayout),
							 (Ā2A1B̄2B1(), GeneralConjugation, GeneralLayout))
		lat = GrassmannLattice(N=2, δτ=0.1, contour=:imag, ordering=o)
		@test lat isa ImagGrassmannLattice
		@test scalartype(lat) == Float64
		@test ConjugationStyle(lat) isa conj
		@test LayoutStyle(lat) isa layout
	end
	for (o, conj, layout) in ((A1Ā1B1B̄1a1ā1b1b̄1(), AdjacentConjugation, TimeLocalLayout),
							  (A1Ā1a1ā1B1B̄1b1b̄1(), AdjacentConjugation, TimeLocalLayout),
							  (A1Ā1B1B̄1b̄1B̄1ā1Ā1(), GeneralConjugation, TimeLocalLayout),
							  (A1B1ā1b̄1Ā1B̄1a1b1(), GeneralConjugation, TimeLocalLayout),
							  (A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1(), AdjacentConjugation, BandLocalLayout),
							  (A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), GeneralConjugation, BranchLocalLayout),
							  (A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2(), AdjacentConjugation, BranchLocalLayout),
							  (Ā2A1ā1a2B̄2B1b̄1b̄2(), GeneralConjugation, GeneralLayout))
		lat = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=o)
		@test lat isa RealGrassmannLattice
		@test scalartype(lat) == ComplexF64
		@test ConjugationStyle(lat) isa conj
		@test LayoutStyle(lat) isa layout
	end
	for (o, conj, layout) in ((A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2(), AdjacentConjugation, TimeLocalLayout),
							  (A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2(), AdjacentConjugation, TimeLocalLayout),
							  (A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), GeneralConjugation, BranchLocalLayout),
							  (Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2(), GeneralConjugation, GeneralLayout),
							  (A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1(), GeneralConjugation, TimeLocalLayout))
		lat = GrassmannLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=o)
		@test lat isa MixedGrassmannLattice
		@test scalartype(lat) == ComplexF64
		@test ConjugationStyle(lat) isa conj
		@test LayoutStyle(lat) isa layout
	end
end

# helper: compare an expected site-position map against index().
# real-time lattices reject branch :τ; their i=0 positions do not depend on the
# branch, so for :τ entries the kwarg is omitted there.
check_indexmap(lat, expected::Dict) = all(
	let p = pos
		(br == :τ && lat isa RealGrassmannLattice) ? index(lat, i, conj=c, band=b) == p :
		index(lat, i, conj=c, branch=br, band=b) == p
	end
	for ((i, c, br, b), pos) in expected)

@testset "API: lattice index maps, imaginary time (all orderings)" begin
	# a\bar{a}b\bar{b} a_3\bar{a}_3b_3\bar{b}_3 a_2\bar{a}_2b_2\bar{b}_2 a_1\bar{a}_1b_1\bar{b}_1
	lat = GrassmannLattice(N=2, δτ=0.05, bands=2, contour=:imag, ordering=A1Ā1B1B̄1())
	@test length(lat) == 16 && lat.k == 3 && lat.β == 0.1 && lat.T == 10
	@test lat.τs == 0:0.05:0.1 && lat.bands == 2 && lat.δτ == 0.05
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(3, false, :τ, 1) => 5,  (3, true, :τ, 1) => 6,  (3, false, :τ, 2) => 7,  (3, true, :τ, 2) => 8,
		(2, false, :τ, 1) => 9,  (2, true, :τ, 1) => 10, (2, false, :τ, 2) => 11, (2, true, :τ, 2) => 12,
		(1, false, :τ, 1) => 13, (1, true, :τ, 1) => 14, (1, false, :τ, 2) => 15, (1, true, :τ, 2) => 16))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# ab\bar{b}\bar{a} a_3b_3\bar{b}_3\bar{a}_3 a_2b_2\bar{b}_2\bar{a}_2 a_1b_1\bar{b}_1\bar{a}_1
	lat = GrassmannLattice(N=2, δτ=0.05, bands=2, contour=:imag, ordering=A1B1B̄1Ā1())
	@test length(lat) == 16 && lat.k == 3 && lat.β == 0.1 && lat.T == 10
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, false, :τ, 2) => 2,  (0, true, :τ, 2) => 3,  (0, true, :τ, 1) => 4,
		(3, false, :τ, 1) => 5,  (3, false, :τ, 2) => 6,  (3, true, :τ, 2) => 7,  (3, true, :τ, 1) => 8,
		(2, false, :τ, 1) => 9,  (2, false, :τ, 2) => 10, (2, true, :τ, 2) => 11, (2, true, :τ, 1) => 12,
		(1, false, :τ, 1) => 13, (1, false, :τ, 2) => 14, (1, true, :τ, 2) => 15, (1, true, :τ, 1) => 16))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# a\bar{a}b\bar{b} a_2\bar{a}_2a_1\bar{a}_1 b_2\bar{b}_2b_1\bar{b}_1
	lat = GrassmannLattice(N=1, δτ=0.05, bands=2, contour=:imag, ordering=A2Ā2A1Ā1B2B̄2B1B̄1())
	@test length(lat) == 12 && lat.k == 2 && lat.β == 0.05 && lat.T == 20
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(2, false, :τ, 1) => 5,  (2, true, :τ, 1) => 6,  (1, false, :τ, 1) => 7,  (1, true, :τ, 1) => 8,
		(2, false, :τ, 2) => 9,  (2, true, :τ, 2) => 10, (1, false, :τ, 2) => 11, (1, true, :τ, 2) => 12))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# aābb̄ a₃b₃ ā₃a₂b̄₃b₂ ā₂a₁b̄₂b₁ ā₁b̄₁ (retarded-interaction ordering)
	lat = GrassmannLattice(N=2, δτ=0.05, bands=2, contour=:imag, ordering=Ā2A1B̄2B1())
	@test length(lat) == 16 && lat.k == 3 && lat.β == 0.1 && lat.T == 10
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(3, false, :τ, 1) => 5,  (3, false, :τ, 2) => 6,
		(3, true, :τ, 1) => 7,   (2, false, :τ, 1) => 8, (3, true, :τ, 2) => 9,   (2, false, :τ, 2) => 10,
		(2, true, :τ, 1) => 11,  (1, false, :τ, 1) => 12, (2, true, :τ, 2) => 13, (1, false, :τ, 2) => 14,
		(1, true, :τ, 1) => 15,  (1, true, :τ, 2) => 16))
end

@testset "API: lattice index maps, real time (all orderings)" begin
	# a\bar{a} a_3^+\bar{a}_3^+a_3^-\bar{a}_3^- a_2^+... a_1^+...
	lat = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A1Ā1B1B̄1a1ā1b1b̄1())
	@test length(lat) == 14 && lat.k == 3 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :+, 1) => 3,  (3, true, :+, 1) => 4,  (3, false, :-, 1) => 5,  (3, true, :-, 1) => 6,
		(2, false, :+, 1) => 7,  (2, true, :+, 1) => 8,  (2, false, :-, 1) => 9,  (2, true, :-, 1) => 10,
		(1, false, :+, 1) => 11, (1, true, :+, 1) => 12, (1, false, :-, 1) => 13, (1, true, :-, 1) => 14))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# aābb̄ a₂+ā₂+a₂-ā₂-b₂+b̄₂+b₂-b̄₂- a₁+...
	lat = GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=A1Ā1a1ā1B1B̄1b1b̄1())
	@test length(lat) == 20 && lat.k == 2 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(2, false, :+, 1) => 5,  (2, true, :+, 1) => 6,  (2, false, :-, 1) => 7,  (2, true, :-, 1) => 8,
		(2, false, :+, 2) => 9,  (2, true, :+, 2) => 10, (2, false, :-, 2) => 11, (2, true, :-, 2) => 12,
		(1, false, :+, 1) => 13, (1, true, :+, 1) => 14, (1, false, :-, 1) => 15, (1, true, :-, 1) => 16,
		(1, false, :+, 2) => 17, (1, true, :+, 2) => 18, (1, false, :-, 2) => 19, (1, true, :-, 2) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# abb̄ā a₂+a₂-b₂+b₂-b̄₂-b̄₂+ā₂-ā₂+ a₁+... (historical ordering)
	lat = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A1Ā1B1B̄1b̄1B̄1ā1Ā1())
	@test length(lat) == 14 && lat.k == 3 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :+, 1) => 3,  (3, false, :-, 1) => 4, (3, true, :-, 1) => 5,  (3, true, :+, 1) => 6,
		(2, false, :+, 1) => 7,  (2, false, :-, 1) => 8, (2, true, :-, 1) => 9,  (2, true, :+, 1) => 10,
		(1, false, :+, 1) => 11, (1, false, :-, 1) => 12, (1, true, :-, 1) => 13, (1, true, :+, 1) => 14))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# a₂+a₂-b₂+b₂-b̄₂-b̄₂+ā₂-ā₂+ with bands
	lat = GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=A1Ā1B1B̄1b̄1B̄1ā1Ā1())
	@test length(lat) == 20 && lat.k == 2 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, false, :τ, 2) => 2,  (0, true, :τ, 2) => 3,  (0, true, :τ, 1) => 4,
		(2, false, :+, 1) => 5,  (2, false, :-, 1) => 6,  (2, false, :+, 2) => 7,  (2, false, :-, 2) => 8,
		(2, true, :-, 2) => 9,   (2, true, :+, 2) => 10,  (2, true, :-, 1) => 11, (2, true, :+, 1) => 12,
		(1, false, :+, 1) => 13, (1, false, :-, 1) => 14, (1, false, :+, 2) => 15, (1, false, :-, 2) => 16,
		(1, true, :-, 2) => 17,  (1, true, :+, 2) => 18,  (1, true, :-, 1) => 19, (1, true, :+, 1) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# aā a₃+ā₃-ā₃+a₃⁻- ... (impurity-dynamics time-local ordering)
	lat = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A1B1ā1b̄1Ā1B̄1a1b1())
	@test length(lat) == 14 && lat.k == 3 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :+, 1) => 3,  (3, true, :-, 1) => 4,  (3, true, :+, 1) => 5,  (3, false, :-, 1) => 6,
		(2, false, :+, 1) => 7,  (2, true, :-, 1) => 8,  (2, true, :+, 1) => 9,  (2, false, :-, 1) => 10,
		(1, false, :+, 1) => 11, (1, true, :-, 1) => 12, (1, true, :+, 1) => 13, (1, false, :-, 1) => 14))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# a₂+b₂+ā₂-b̄₂-ā₂+b̄₂+a₂-b₂- with bands
	lat = GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=A1B1ā1b̄1Ā1B̄1a1b1())
	@test length(lat) == 20 && lat.k == 2 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(2, false, :+, 1) => 5,  (2, false, :+, 2) => 6,  (2, true, :-, 1) => 7,  (2, true, :-, 2) => 8,
		(2, true, :+, 1) => 9,   (2, true, :+, 2) => 10,  (2, false, :-, 1) => 11, (2, false, :-, 2) => 12,
		(1, false, :+, 1) => 13, (1, false, :+, 2) => 14, (1, true, :-, 1) => 15, (1, true, :-, 2) => 16,
		(1, true, :+, 1) => 17,  (1, true, :+, 2) => 18,  (1, false, :-, 1) => 19, (1, false, :-, 2) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# aābb̄ a₂+ā₂+a₁+ā₁+ a₂-ā₂-a₁-ā₁- b₂+b̄₂+b₁+b̄₁+ b₂-b̄₂-b₁-b̄₁-
	lat = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1())
	@test length(lat) == 14 && lat.k == 3 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :+, 1) => 3,  (3, true, :+, 1) => 4,  (2, false, :+, 1) => 5,  (2, true, :+, 1) => 6,
		(1, false, :+, 1) => 7,  (1, true, :+, 1) => 8,
		(3, false, :-, 1) => 9,  (3, true, :-, 1) => 10, (2, false, :-, 1) => 11, (2, true, :-, 1) => 12,
		(1, false, :-, 1) => 13, (1, true, :-, 1) => 14))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# band-local ordering with bands
	lat = GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1())
	@test length(lat) == 20 && lat.k == 2 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(2, false, :+, 1) => 5,  (2, true, :+, 1) => 6,  (1, false, :+, 1) => 7,  (1, true, :+, 1) => 8,
		(2, false, :-, 1) => 9,  (2, true, :-, 1) => 10, (1, false, :-, 1) => 11, (1, true, :-, 1) => 12,
		(2, false, :+, 2) => 13, (2, true, :+, 2) => 14, (1, false, :+, 2) => 15, (1, true, :+, 2) => 16,
		(2, false, :-, 2) => 17, (2, true, :-, 2) => 18, (1, false, :-, 2) => 19, (1, true, :-, 2) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# abb̄ā a₂+b₂+b̄₂+ā₂+a₁+b₁+b̄₁+ā₁+ a₁-b₁-b̄₁-ā₁-a₂-b₂-b̄₂-ā₂-
	lat = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2())
	@test length(lat) == 14 && lat.k == 3 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :+, 1) => 3,  (3, true, :+, 1) => 4,  (2, false, :+, 1) => 5,  (2, true, :+, 1) => 6,
		(1, false, :+, 1) => 7,  (1, true, :+, 1) => 8,
		(1, false, :-, 1) => 9,  (1, true, :-, 1) => 10, (2, false, :-, 1) => 11, (2, true, :-, 1) => 12,
		(3, false, :-, 1) => 13, (3, true, :-, 1) => 14))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# branch-local ordering with bands
	lat = GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2())
	@test length(lat) == 20 && lat.k == 2 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, false, :τ, 2) => 2,  (0, true, :τ, 2) => 3,  (0, true, :τ, 1) => 4,
		(2, false, :+, 1) => 5,  (2, false, :+, 2) => 6,  (2, true, :+, 2) => 7,  (2, true, :+, 1) => 8,
		(1, false, :+, 1) => 9,  (1, false, :+, 2) => 10, (1, true, :+, 2) => 11, (1, true, :+, 1) => 12,
		(1, false, :-, 1) => 13, (1, false, :-, 2) => 14, (1, true, :-, 2) => 15, (1, true, :-, 1) => 16,
		(2, false, :-, 1) => 17, (2, false, :-, 2) => 18, (2, true, :-, 2) => 19, (2, true, :-, 1) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# aābb̄ a₂+ā₂+b₂+b̄₂+a₁+ā₁+b₁+b̄₁+ a₁-ā₁-b₁-b̄₁-a₂-ā₂-b₂-b̄₂-
	lat = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
	@test length(lat) == 14 && lat.k == 3 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :+, 1) => 3,  (3, true, :+, 1) => 4,  (2, false, :+, 1) => 5,  (2, true, :+, 1) => 6,
		(1, false, :+, 1) => 7,  (1, true, :+, 1) => 8,
		(1, false, :-, 1) => 9,  (1, true, :-, 1) => 10, (2, false, :-, 1) => 11, (2, true, :-, 1) => 12,
		(3, false, :-, 1) => 13, (3, true, :-, 1) => 14))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# branch-local ordering with bands
	lat = GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
	@test length(lat) == 20 && lat.k == 2 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(2, false, :+, 1) => 5,  (2, true, :+, 1) => 6,  (2, false, :+, 2) => 7,  (2, true, :+, 2) => 8,
		(1, false, :+, 1) => 9,  (1, true, :+, 1) => 10, (1, false, :+, 2) => 11, (1, true, :+, 2) => 12,
		(1, false, :-, 1) => 13, (1, true, :-, 1) => 14, (1, false, :-, 2) => 15, (1, true, :-, 2) => 16,
		(2, false, :-, 1) => 17, (2, true, :-, 1) => 18, (2, false, :-, 2) => 19, (2, true, :-, 2) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# aābb̄ a₃+ā₃-b₃+b̄₃- ā₃+a₂+ā₂-a₃-b̄₃+b₂+b̄₂-b₃- ... (retarded interaction)
	lat = GrassmannLattice(N=2, δt=0.05, contour=:real, ordering=Ā2A1ā1a2B̄2B1b̄1b̄2())
	@test length(lat) == 14 && lat.k == 3 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :+, 1) => 3,  (3, true, :-, 1) => 4,
		(3, true, :+, 1) => 5,   (2, false, :+, 1) => 6, (2, true, :-, 1) => 7,  (3, false, :-, 1) => 8,
		(2, true, :+, 1) => 9,   (1, false, :+, 1) => 10, (1, true, :-, 1) => 11, (2, false, :-, 1) => 12,
		(1, true, :+, 1) => 13,  (1, false, :-, 1) => 14))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# retarded-interaction ordering with bands
	lat = GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=Ā2A1ā1a2B̄2B1b̄1b̄2())
	@test length(lat) == 20 && lat.k == 2 && lat.t == 0.1
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(2, false, :+, 1) => 5,  (2, true, :-, 1) => 6,  (2, false, :+, 2) => 7,  (2, true, :-, 2) => 8,
		(2, true, :+, 1) => 9,   (1, false, :+, 1) => 10, (1, true, :-, 1) => 11, (2, false, :-, 1) => 12,
		(2, true, :+, 2) => 13,  (1, false, :+, 2) => 14, (1, true, :-, 2) => 15, (2, false, :-, 2) => 16,
		(1, true, :+, 1) => 17,  (1, false, :-, 1) => 18, (1, true, :+, 2) => 19, (1, false, :-, 2) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6
end

@testset "API: lattice index maps, mixed time (all orderings)" begin
	# imag: aābb̄ a₂ā₂b₂b̄₂ a₁ā₁b₁b̄₁ / real: a₃+ā₃+a₃-ā₃- ... ascending in time
	lat = GrassmannLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed,
						   ordering=A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2())
	@test length(lat) == 20 && lat.Nt == 2 && lat.Nτ == 2
	@test lat.t == 0.1 && lat.β == 0.2 && lat.ts == 0:0.05:0.1 && lat.τs == 0:0.1:0.2
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :τ, 1) => 3,  (3, true, :τ, 1) => 4,  (2, false, :τ, 1) => 5,  (2, true, :τ, 1) => 6,
		(1, false, :τ, 1) => 7,  (1, true, :τ, 1) => 8,
		(1, false, :+, 1) => 9,  (1, true, :+, 1) => 10, (1, false, :-, 1) => 11, (1, true, :-, 1) => 12,
		(2, false, :+, 1) => 13, (2, true, :+, 1) => 14, (2, false, :-, 1) => 15, (2, true, :-, 1) => 16,
		(3, false, :+, 1) => 17, (3, true, :+, 1) => 18, (3, false, :-, 1) => 19, (3, true, :-, 1) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	lat = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, bands=2, contour=:mixed,
						   ordering=A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2())
	@test length(lat) == 32 && lat.Nt == 1 && lat.Nτ == 2
	@test lat.t == 0.1 && lat.β == 0.2
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(3, false, :τ, 1) => 5,  (3, true, :τ, 1) => 6,  (3, false, :τ, 2) => 7,  (3, true, :τ, 2) => 8,
		(2, false, :τ, 1) => 9,  (2, true, :τ, 1) => 10, (2, false, :τ, 2) => 11, (2, true, :τ, 2) => 12,
		(1, false, :τ, 1) => 13, (1, true, :τ, 1) => 14, (1, false, :τ, 2) => 15, (1, true, :τ, 2) => 16,
		(1, false, :+, 1) => 17, (1, true, :+, 1) => 18, (1, false, :-, 1) => 19, (1, true, :-, 1) => 20,
		(1, false, :+, 2) => 21, (1, true, :+, 2) => 22, (1, false, :-, 2) => 23, (1, true, :-, 2) => 24,
		(2, false, :+, 1) => 25, (2, true, :+, 1) => 26, (2, false, :-, 1) => 27, (2, true, :-, 1) => 28,
		(2, false, :+, 2) => 29, (2, true, :+, 2) => 30, (2, false, :-, 2) => 31, (2, true, :-, 2) => 32))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# imag: aābb̄ / real: a₃-ā₃-... ascending (branches swapped)
	lat = GrassmannLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed,
						   ordering=A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2())
	@test length(lat) == 20 && lat.Nt == 2 && lat.Nτ == 2
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :τ, 1) => 3,  (3, true, :τ, 1) => 4,  (2, false, :τ, 1) => 5,  (2, true, :τ, 1) => 6,
		(1, false, :τ, 1) => 7,  (1, true, :τ, 1) => 8,
		(1, false, :-, 1) => 9,  (1, true, :-, 1) => 10, (1, false, :+, 1) => 11, (1, true, :+, 1) => 12,
		(2, false, :-, 1) => 13, (2, true, :-, 1) => 14, (2, false, :+, 1) => 15, (2, true, :+, 1) => 16,
		(3, false, :-, 1) => 17, (3, true, :-, 1) => 18, (3, false, :+, 1) => 19, (3, true, :+, 1) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	lat = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, bands=2, contour=:mixed,
						   ordering=A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2())
	@test length(lat) == 32 && lat.Nt == 1 && lat.Nτ == 2
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(3, false, :τ, 1) => 5,  (3, true, :τ, 1) => 6,  (3, false, :τ, 2) => 7,  (3, true, :τ, 2) => 8,
		(2, false, :τ, 1) => 9,  (2, true, :τ, 1) => 10, (2, false, :τ, 2) => 11, (2, true, :τ, 2) => 12,
		(1, false, :τ, 1) => 13, (1, true, :τ, 1) => 14, (1, false, :τ, 2) => 15, (1, true, :τ, 2) => 16,
		(1, false, :-, 1) => 17, (1, true, :-, 1) => 18, (1, false, :+, 1) => 19, (1, true, :+, 1) => 20,
		(1, false, :-, 2) => 21, (1, true, :-, 2) => 22, (1, false, :+, 2) => 23, (1, true, :+, 2) => 24,
		(2, false, :-, 1) => 25, (2, true, :-, 1) => 26, (2, false, :+, 1) => 27, (2, true, :+, 1) => 28,
		(2, false, :-, 2) => 29, (2, true, :-, 2) => 30, (2, false, :+, 2) => 31, (2, true, :+, 2) => 32))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# imag: abb̄ā / real: a₂+b₂+b̄₂+ā₂+... (branch-local)
	lat = GrassmannLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed,
						   ordering=A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2())
	@test length(lat) == 20 && lat.Nt == 2 && lat.Nτ == 2
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :τ, 1) => 3,  (3, true, :τ, 1) => 4,  (2, false, :τ, 1) => 5,  (2, true, :τ, 1) => 6,
		(1, false, :τ, 1) => 7,  (1, true, :τ, 1) => 8,
		(3, false, :+, 1) => 9,  (3, true, :+, 1) => 10, (2, false, :+, 1) => 11, (2, true, :+, 1) => 12,
		(1, false, :+, 1) => 13, (1, true, :+, 1) => 14,
		(1, false, :-, 1) => 15, (1, true, :-, 1) => 16, (2, false, :-, 1) => 17, (2, true, :-, 1) => 18,
		(3, false, :-, 1) => 19, (3, true, :-, 1) => 20))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	lat = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, bands=2, contour=:mixed,
						   ordering=A1B1B̄1Ā1_A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2())
	@test length(lat) == 32 && lat.Nt == 1 && lat.Nτ == 2
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, false, :τ, 2) => 2,  (0, true, :τ, 2) => 3,  (0, true, :τ, 1) => 4,
		(3, false, :τ, 1) => 5,  (3, false, :τ, 2) => 6,  (3, true, :τ, 2) => 7,  (3, true, :τ, 1) => 8,
		(2, false, :τ, 1) => 9,  (2, false, :τ, 2) => 10, (2, true, :τ, 2) => 11, (2, true, :τ, 1) => 12,
		(1, false, :τ, 1) => 13, (1, false, :τ, 2) => 14, (1, true, :τ, 2) => 15, (1, true, :τ, 1) => 16,
		(2, false, :+, 1) => 17, (2, false, :+, 2) => 18, (2, true, :+, 2) => 19, (2, true, :+, 1) => 20,
		(1, false, :+, 1) => 21, (1, false, :+, 2) => 22, (1, true, :+, 2) => 23, (1, true, :+, 1) => 24,
		(1, false, :-, 1) => 25, (1, false, :-, 2) => 26, (1, true, :-, 2) => 27, (1, true, :-, 1) => 28,
		(2, false, :-, 1) => 29, (2, false, :-, 2) => 30, (2, true, :-, 2) => 31, (2, true, :-, 1) => 32))
	@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6

	# retarded-interaction mixed-time ordering
	lat = GrassmannLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed,
						   ordering=Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2())
	@test length(lat) == 20 && lat.Nt == 2 && lat.Nτ == 2
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,
		(3, false, :τ, 1) => 3,  (1, true, :τ, 1) => 4,
		(1, false, :-, 1) => 5,  (1, true, :+, 1) => 6,
		(3, true, :τ, 1) => 7,   (2, false, :τ, 1) => 8, (2, true, :τ, 1) => 9,  (1, false, :τ, 1) => 10,
		(1, true, :-, 1) => 11,  (2, false, :-, 1) => 12, (2, true, :+, 1) => 13, (1, false, :+, 1) => 14,
		(2, true, :-, 1) => 15,  (3, false, :-, 1) => 16, (3, true, :+, 1) => 17, (2, false, :+, 1) => 18,
		(3, true, :-, 1) => 19,  (3, false, :+, 1) => 20))

	lat = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, bands=2, contour=:mixed,
						   ordering=Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2())
	@test length(lat) == 32 && lat.Nt == 1 && lat.Nτ == 2
	@test check_indexmap(lat, Dict(
		(0, false, :τ, 1) => 1,  (0, true, :τ, 1) => 2,  (0, false, :τ, 2) => 3,  (0, true, :τ, 2) => 4,
		(3, false, :τ, 1) => 5,  (3, false, :τ, 2) => 6,
		(1, true, :τ, 1) => 7,   (1, true, :τ, 2) => 8,
		(1, false, :-, 1) => 9,  (1, true, :+, 1) => 10, (1, false, :-, 2) => 11, (1, true, :+, 2) => 12,
		(3, true, :τ, 1) => 13,  (2, false, :τ, 1) => 14, (3, true, :τ, 2) => 15, (2, false, :τ, 2) => 16,
		(2, true, :τ, 1) => 17,  (1, false, :τ, 1) => 18, (2, true, :τ, 2) => 19, (1, false, :τ, 2) => 20,
		(1, true, :-, 1) => 21,  (2, false, :-, 1) => 22, (2, true, :+, 1) => 23, (1, false, :+, 1) => 24,
		(1, true, :-, 2) => 25,  (2, false, :-, 2) => 26, (2, true, :+, 2) => 27, (1, false, :+, 2) => 28,
		(2, true, :-, 1) => 29,  (2, false, :+, 1) => 30, (2, true, :-, 2) => 31, (2, false, :+, 2) => 32))
end

# indexmappings must agree with index() and be a bijection for every ordering
@testset "API: indexmappings consistency (all orderings)" begin
	for o in all_imag_orderings
		lat = GrassmannLattice(N=2, δτ=0.1, bands=2, contour=:imag, ordering=o)
		m = indexmappings(lat)
		@test length(m) == length(lat) == 2 * lat.bands * (lat.k + 1)
		@test length(unique(values(m))) == length(m)
		@test all(m[(i, c, :τ, b)] == index(lat, i, conj=c, band=b)
				  for i in 0:lat.k for c in (true, false) for b in 1:lat.bands)
	end
	for o in all_real_orderings
		lat = GrassmannLattice(N=2, δt=0.1, bands=2, contour=:real, ordering=o)
		m = indexmappings(lat)
		@test length(m) == length(lat) == 4 * lat.bands * lat.k + 2 * lat.bands
		@test length(unique(values(m))) == length(m)
		@test all(m[(0, c, :+, b)] == index(lat, 0, conj=c, branch=:+, band=b)
				  for c in (true, false) for b in 1:lat.bands)
		@test all(m[(i, c, br, b)] == index(lat, i, conj=c, branch=br, band=b)
				  for i in 1:lat.k for c in (true, false) for br in (:+, :-) for b in 1:lat.bands)
	end
	for o in all_mixed_orderings
		lat = GrassmannLattice(Nt=2, δt=0.1, Nτ=2, δτ=0.1, bands=2, contour=:mixed, ordering=o)
		m = indexmappings(lat)
		@test length(m) == length(lat)
		@test length(unique(values(m))) == length(m)
		@test all(m[(i, c, :τ, b)] == index(lat, i, conj=c, branch=:τ, band=b)
				  for i in 1:lat.kτ for c in (true, false) for b in 1:lat.bands)
		@test all(m[(0, c, :+, b)] == index(lat, 0, conj=c, branch=:+, band=b)
				  for c in (true, false) for b in 1:lat.bands)
		@test all(m[(i, c, br, b)] == index(lat, i, conj=c, branch=br, band=b)
				  for i in 1:lat.kt for c in (true, false) for br in (:+, :-) for b in 1:lat.bands)
	end
end

# integration of the uniform and the vacuum states over the whole lattice gives 1.
# skipped for the two GeneralConjugation mixed-time orderings: toadjacentordering
# is not implemented for them (consistent with the archived tests).
@testset "API: integrate uniform & vacuum states (all orderings)" begin
	for o in all_imag_orderings
		lat = GrassmannLattice(N=2, δτ=0.05, bands=2, contour=:imag, ordering=o)
		@test integrate(lat, GrassmannMPS(scalartype(lat), length(lat))) ≈ 1 atol = 1.0e-6
		@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6
	end
	for o in all_real_orderings
		lat = GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=o)
		@test integrate(lat, GrassmannMPS(scalartype(lat), length(lat))) ≈ 1 atol = 1.0e-6
		@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6
	end
	for o in all_mixed_orderings
		(typeof(o) <: Union{typeof(Ā3A2B̄3B2Ā2A1B̄2B1_ā1a2Ā2A1b̄1b2B̄2B1ā2a3Ā3A2b̄2b3B̄3B2()),
							typeof(A1B1B̄1Ā1_a1b1Ā1B̄1ā1b̄1A1B1())}) && continue
		lat = GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, bands=2, contour=:mixed, ordering=o)
		@test integrate(lat, GrassmannMPS(scalartype(lat), length(lat))) ≈ 1 atol = 1.0e-6
		@test integrate(lat, vacuumstate(lat)) ≈ 1 atol = 1.0e-6
	end
end

# changeordering: matchindices must give a bijection between every pair of orderings
@testset "API: changeordering matchindices (all ordering pairs)" begin
	for (orderings, mk) in ((all_imag_orderings, o -> GrassmannLattice(N=2, δτ=0.1, bands=2, contour=:imag, ordering=o)),
							(all_real_orderings, o -> GrassmannLattice(N=1, δt=0.1, bands=2, contour=:real, ordering=o)),
							(all_mixed_orderings, o -> GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.03, bands=2, contour=:mixed, ordering=o)))
		for o1 in orderings, o2 in orderings
			tsc, src = mk(o1), mk(o2)
			mapping = matchindices(tsc, src)
			@test sort!([mapping[i] for i in 1:length(mapping)]) == collect(1:length(mapping))
		end
	end
end

# swapbandperm must be consistent with indexmappings for every ordering
@testset "API: swapbandperm (all orderings)" begin
	for (orderings, mk) in ((all_imag_orderings, o -> GrassmannLattice(N=2, δτ=0.1, bands=3, contour=:imag, ordering=o)),
							(all_real_orderings, o -> GrassmannLattice(N=2, δt=0.1, bands=3, contour=:real, ordering=o)),
							(all_mixed_orderings, o -> GrassmannLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.03, bands=3, contour=:mixed, ordering=o)))
		for o in orderings
			x = mk(o)
			m = indexmappings(x)
			for (b1, b2) in ((1, 2), (2, 3))
				perm = swapbandperm(x, b1, b2)
				change_band(b) = ifelse(b == b1, b2, ifelse(b == b2, b1, b))
				@test all(m[(j, c, br, change_band(band))] == perm[pos] for ((j, c, br, band), pos) in m)
			end
		end
	end
end

@testset "API: FockLattice (electron-phonon path)" begin
	# imaginary time M1N1
	lat = FockLattice(N=1, δτ=0.1, contour=:imag, ordering=M1N1())
	@test lat isa ImagFockLattice
	@test LayoutStyle(lat) isa TimeLocalLayout
	@test scalartype(lat) == Float64
	@test length(lat) == 1 && lat.N == 1 && lat.β == 0.1 && lat.T == 10 && lat.bands == 1
	@test lat.τs == 0:0.1:0.1
	@test index(lat, 1) == 1

	lat = FockLattice(N=2, δτ=0.05, bands=2, contour=:imag, ordering=M1N1())
	@test lat isa ImagFockLattice
	@test length(lat) == 4 && lat.N == 2 && lat.β == 0.1 && lat.bands == 2 && lat.T == 10
	@test index(lat, 2, band=1) == 1 && index(lat, 2, band=2) == 2
	@test index(lat, 1, band=1) == 3 && index(lat, 1, band=2) == 4

	lat = FockLattice(N=2, δτ=0.1, bands=3, contour=:imag, ordering=M1N1())
	@test length(lat) == 6 && lat.N == 2 && lat.β == 0.2 && lat.bands == 3 && lat.τs == 0:0.1:0.2
	@test all(index(lat, 2, band=b) == b && index(lat, 1, band=b) == 3 + b for b in 1:3)

	# real time M1m1N1n1
	lat = FockLattice(N=1, δt=0.1, contour=:real, ordering=M1m1N1n1())
	@test lat isa RealFockLattice
	@test LayoutStyle(lat) isa TimeLocalLayout
	@test scalartype(lat) == ComplexF64
	@test length(lat) == 2 && lat.N == 1 && lat.t == 0.1 && lat.δt == 0.1
	@test index(lat, 1, branch=:+) == 1 && index(lat, 1, branch=:-) == 2

	lat = FockLattice(N=2, δt=0.05, bands=2, contour=:real, ordering=M1m1N1n1())
	@test length(lat) == 8 && lat.N == 2 && lat.t == 0.1 && lat.bands == 2
	@test index(lat, 2, band=1, branch=:+) == 1 && index(lat, 2, band=1, branch=:-) == 2
	@test index(lat, 2, band=2, branch=:+) == 3 && index(lat, 2, band=2, branch=:-) == 4
	@test index(lat, 1, band=1, branch=:+) == 5 && index(lat, 1, band=1, branch=:-) == 6
	@test index(lat, 1, band=2, branch=:+) == 7 && index(lat, 1, band=2, branch=:-) == 8

	lat = FockLattice(N=2, δt=0.1, bands=3, contour=:real, ordering=M1m1N1n1())
	@test length(lat) == 12 && lat.N == 2 && lat.t == 0.2 && lat.bands == 3 && lat.ts == 0:0.1:0.2
	@test all(index(lat, 2, band=b, branch=:+) == 2*(b-1) + 1 && index(lat, 2, band=b, branch=:-) == 2*b for b in 1:3)
	@test all(index(lat, 1, band=b, branch=:+) == 6 + 2*(b-1) + 1 && index(lat, 1, band=b, branch=:-) == 6 + 2*b for b in 1:3)

	# mixed time M1N1_m1M1n1N1m2M2n2N2
	lat = FockLattice(Nt=2, δt=0.05, Nτ=2, δτ=0.1, contour=:mixed, ordering=M1N1_m1M1n1N1m2M2n2N2())
	@test lat isa MixedFockLattice
	@test LayoutStyle(lat) isa TimeLocalLayout
	@test scalartype(lat) == ComplexF64
	@test length(lat) == 6 && lat.Nt == 2 && lat.Nτ == 2
	@test lat.t == 0.1 && lat.β == 0.2 && lat.ts == 0:0.05:0.1 && lat.τs == 0:0.1:0.2
	@test index(lat, 2, branch=:τ) == 1 && index(lat, 1, branch=:τ) == 2
	@test index(lat, 1, branch=:-) == 3 && index(lat, 1, branch=:+) == 4
	@test index(lat, 2, branch=:-) == 5 && index(lat, 2, branch=:+) == 6

	lat = FockLattice(Nt=1, δt=0.1, Nτ=2, δτ=0.1, bands=2, contour=:mixed, ordering=M1N1_m1M1n1N1m2M2n2N2())
	@test length(lat) == 8 && lat.Nt == 1 && lat.Nτ == 2 && lat.t == 0.1 && lat.β == 0.2
	@test index(lat, 2, branch=:τ, band=1) == 1 && index(lat, 2, branch=:τ, band=2) == 2
	@test index(lat, 1, branch=:τ, band=1) == 3 && index(lat, 1, branch=:τ, band=2) == 4
	@test index(lat, 1, branch=:-, band=1) == 5 && index(lat, 1, branch=:+, band=1) == 6
	@test index(lat, 1, branch=:-, band=2) == 7 && index(lat, 1, branch=:+, band=2) == 8

	lat = FockLattice(Nt=1, δt=0.1, Nτ=1, δτ=0.1, bands=3, contour=:mixed, ordering=M1N1_m1M1n1N1m2M2n2N2())
	@test length(lat) == 9 && lat.Nt == 1 && lat.Nτ == 1 && lat.t == 0.1 && lat.β == 0.1
	@test all(index(lat, 1, branch=:τ, band=b) == b for b in 1:3)
	@test all(index(lat, 1, branch=:-, band=b) == 2*b + 2 && index(lat, 1, branch=:+, band=b) == 2*b + 3 for b in 1:3)
end
