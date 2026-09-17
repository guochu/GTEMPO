using LinearAlgebra: normalize!

@testset "API: GrassmannMPS" begin
	L, D = 6, 8

	# --- construction and accessors ---
	for T in (Float64, ComplexF64)
		x = randomgmps(T, L, D=D)
		@test scalartype(x) == T
		@test spacetype(x) == Z2Space
		@test length(x) == L
		@test !isempty(x)
		@test firstindex(x) == 1 && lastindex(x) == L
		@test bond_dimensions(x) isa Vector{Int} && length(bond_dimensions(x)) == L
		@test bond_dimension(x) == maximum(bond_dimensions(x))
		@test bond_dimensions(x)[end] == 1
		# virtual spaces chain: space_l(x[i+1]) == space_r(x[i])'
		for i in 1:L-1
			@test space_l(x[i+1]) == space_r(x[i])'
		end
		@test space_l(x[1]) == Z2Space(0 => 1) && space_r(x[L])' == Z2Space(0 => 1)
	end

	# --- svectors bookkeeping ---
	x = randomgmps(L, D=D)
	@test svectors_uninitialized(x)
	@test !iscanonical(x)  # uninitialized svectors never count as mixed-canonical
	unset_svectors!(x)
	@test svectors_uninitialized(x)

	# --- copy / deepcopy independence ---
	y = deepcopy(x)
	y[1] = y[1] * 2
	@test distance(y, x) > 0
	y = copy(x)
	@test distance(y, x) == 0

	# --- complex promotion ---
	xc = complex(x)
	@test scalartype(xc) == ComplexF64
	@test distance(xc, complex(x)) == 0
	@test norm(xc) ≈ norm(x)

	# --- scaling bookkeeping ---
	λ0 = scaling(x)
	nx = norm(x)
	setscaling!(x, 2λ0)
	@test scaling(x) == 2λ0
	@test norm(x) ≈ 2^L * nx  # norm carries scaling^L
	setscaling!(x, λ0)
	@test norm(x) ≈ nx

	# --- norm / dot / distance algebra ---
	x = randomgmps(L, D=D)
	y = randomgmps(L, D=D)
	nx = norm(x)
	ny = norm(y)
	@test norm(x * 2.0) ≈ 2nx
	@test norm(0.5 * x) ≈ 0.5nx
	# complex scalar multiplication on a complex-typed state (a real-typed state
	# cannot absorb a complex factor in place)
	xc = complex(x)
	@test norm(xc * (1.0 + 2.0im)) ≈ abs(1.0 + 2.0im) * norm(xc)
	@test norm(-x) ≈ nx
	@test distance(x / 3.0, x * (1 / 3.0)) / nx < 1.0e-12
	@test dot(x, x) ≈ norm(x)^2
	@test abs(dot(x, y)) <= norm(x) * norm(y) + 1.0e-10  # Cauchy–Schwarz
	@test distance(x, x) == 0
	@test distance(x, 2x) ≈ nx  # |2x - x| = |x|
	lmul!(3.0, x)
	@test norm(x) ≈ 3nx

	# --- linear superposition ---
	x = randomgmps(L, D=D)
	y = randomgmps(L, D=D)
	nx = norm(x)
	@test norm(x + x) ≈ 2 * norm(x) rtol = 1.0e-8
	@test norm(x - x) < 1.0e-8 * nx
	@test distance(x + y, y + x) / nx < 1.0e-8

	# --- canonicalize! / normalize! and canonicality predicates ---
	x = randomgmps(L, D=D)
	canonicalize!(x, alg=Orthogonalize(trunc=truncdimcutoff(D=64, ϵ=1.0e-14), normalize=false))
	@test !svectors_uninitialized(x)
	@test isrightcanonical(x)
	@test iscanonical(x)
	normalize!(x)
	@test norm(x) ≈ 1.0 rtol = 1.0e-10

	# --- increase_bond! ---
	x = randomgmps(L, D=4)
	x0 = deepcopy(x)
	xd = increase_bond!(x, 16)
	@test minimum(bond_dimensions(xd)[2:end-1]) >= 4  # original internal bonds have dim 4
	@test maximum(bond_dimensions(xd)) == 16
	@test svectors_uninitialized(xd)
	@test distance(xd, x0) / norm(x0) < 1.0e-12  # isometry insertion preserves the state

	# --- permute!: stays mixed-canonical, inverse permute restores the state ---
	truncbig = truncdimcutoff(D=256, ϵ=1.0e-12)
	for Lp in (6, 8), Dp in (4, 8)
		x = randomgmps(Lp, D=Dp)
		canonicalize!(x, alg=Orthogonalize(trunc=truncbig, normalize=false))
		@test iscanonical(x)
		nxp = norm(x)

		perms = (
			collect(1:Lp),                        # identity
			[2; 1; collect(3:Lp)],                # single adjacent swap
			reverse(collect(1:Lp)),               # full reversal
			randperm(Lp),                         # random permutation
		)


	# a mixed-canonical random state (right-canonical with Schmidt svectors)
	canonicalized = randomgmps(4, D=4)
	canonicalize!(canonicalized, alg=Orthogonalize(trunc=truncdimcutoff(D=64, ϵ=1.0e-14), normalize=false))

	# --- GTerm / ExpGTerm / z2space / SparseGMPS / togmps / sparse mult! ---
	truncbig = truncdimcutoff(D=256, ϵ=1.0e-12)
	@test z2space() == Z2Space(0 => 1, 1 => 1)

	g = GTerm(3, 1; coeff=2.0)
	@test g isa GTerm && g isa AbstractGTerm
	@test GTEMPO.positions(g) == (1, 3)   # positions get sorted on construction
	@test g.coeff == -2.0          # fermionic sign of the sorting swap
	@test_throws ArgumentError GTerm(1; coeff=1.0)      # odd number of variables
	@test_throws ArgumentError GTerm(2, 2; coeff=1.0)   # duplicated positions

	s = SparseGMPS(GTerm(1, 3; coeff=2.0))
	@test GTEMPO.positions(s) == [1, 2, 3]
	gm = togmps(s, 6)
	@test gm isa GrassmannMPS && length(gm) == 6
	@test norm(gm) ≈ 2.0 atol = 1.0e-12
	@test !isleftcanonical(canonicalized)  # a mixed-canonical random state is not left-canonical

	# in-place mult! with a SparseGMPS agrees with the dense version
	x6 = randomgmps(L, D=4)
	xa = deepcopy(x6); xb = deepcopy(x6)
	mult!(xa, s, trunc=truncbig)
	mult!(xb, togmps(s, 6), trunc=truncbig)
	@test distance(xa, xb) / norm(xa) < 1.0e-6

	# ExpGTerm
	eg = exp(GTerm(1, 3; coeff=0.5))
	@test eg isa ExpGTerm && GTEMPO.positions(eg) == (1, 3)
	@test convert(PartialMPO, eg) isa PartialMPO

	# PartialMPO conversion and in-place apply!
	pmpo = convert(PartialMPO, GTerm(2, 4; coeff=1.0))
	@test pmpo isa PartialMPO
	@test GTEMPO.positions(pmpo) == [2, 3, 4]
	@test distance(pmpo * deepcopy(x), apply!(pmpo, deepcopy(x))) == 0

	# GrassmannTransferMatrix
	tm = GrassmannTransferMatrix(x6, x6)
	@test length(tm) == L
	@test scaling(tm) ≈ (scaling(x6) * scaling(x6))^2
	@test length(GrassmannTransferMatrix(2, x6, x6)) == 2

	# tensor type aliases and abstract types
	@test x[1] isa MPSTensor
	@test x isa AbstractFiniteGMPS && x isa AbstractGMPS
	m4 = zeros(Float64, z2space() ⊗ z2space() ← z2space() ⊗ z2space())
	@test m4 isa MPOTensor && m4 isa SiteOperator
	@test physical_space(m4) == z2space()
	@test iphysical_space(m4) == z2space()'
	@test ophysical_space(m4) == z2space()
	@test physical_space(x, 3) == physical_space(x[3])
	@test physical_spaces(x) == [physical_space(x[i]) for i in 1:length(x)]

	# NoTruncation keeps every Schmidt value
	xn = randomgmps(4, D=4)
	canonicalize!(xn, alg=Orthogonalize(trunc=NoTruncation(), normalize=false))
	@test iscanonical(xn) && bond_dimension(xn) == 4

		for p in perms
			# non-mutating permute: result stays mixed-canonical
			y = permute(x, p; trunc=truncbig)
			@test iscanonical(y)

			# in-place permute agrees with the non-mutating version
			x2 = deepcopy(x)
			permute!(x2, p; trunc=truncbig)
			@test iscanonical(x2)
			@test distance(x2, y) / nxp < 1.0e-6

			# inverse permute restores the original MPS
			w = permute(y, invperm(p); trunc=truncbig)
			@test iscanonical(w)
			@test distance(w, x) / nxp < 1.0e-6
		end
	end
end
