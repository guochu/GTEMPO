# Tests for the @grassmann macro: it shares the @tensor syntax/semantics but
# dispatches every contraction / addition / trace to the fermionic
# (Grassmann) backend. On even-parity tensors it must reproduce the bosonic
# @tensor results exactly; differences only arise from closed fermionic
# (odd-parity) U-turn loops, which contribute a factor -1 each.

@testset "API: @grassmann vs @tensor (even parity)" begin
	s0 = Z2Space(0 => 1)

	# plain contraction
	a = randn(Float64, s0 ⊗ s0 ← s0)
	b = randn(Float64, s0 ← s0 ⊗ s0)
	c1 = zeros(Float64, s0 ⊗ s0 ← s0 ⊗ s0)
	c2 = zeros(Float64, s0 ⊗ s0 ← s0 ⊗ s0)
	@grassmann c1[1, 2; 4, 5] := a[1, 2, 3] * b[3, 4, 5]
	@tensor c2[1, 2; 4, 5] := a[1, 2, 3] * b[3, 4, 5]
	@test norm(c1 - c2) == 0

	# scalar (fully contracted) trace
	A = randn(Float64, s0 ⊗ s0 ← s0 ⊗ s0)
	@test (@grassmann A[1, 2, 1, 2]) ≈ (@tensor A[1, 2, 1, 2]) atol=0.0

	# in-place addition
	d1 = randn(Float64, s0 ⊗ s0 ← s0 ⊗ s0)
	d2 = copy(d1)
	@grassmann d1[1, 2; 4, 5] += a[1, 2, 3] * b[3, 4, 5]
	@tensor d2[1, 2; 4, 5] += a[1, 2, 3] * b[3, 4, 5]
	@test norm(d1 - d2) == 0

	# three-tensor chain with intermediate rank (1, 1) tensors
	t1 = isomorphism(Float64, s0 ⊗ s0, s0)
	t2 = isomorphism(Float64, s0, s0 ⊗ s0)
	e1 = zeros(Float64, s0 ⊗ s0 ← s0)
	e2 = zeros(Float64, s0 ⊗ s0 ← s0)
	@grassmann e1[1, 2; 3] := a[1, 2, 4] * t2[4, 5, 6] * t1[6, 5, 3]
	@tensor e2[1, 2; 3] := a[1, 2, 4] * t2[4, 5, 6] * t1[6, 5, 3]
	@test norm(e1 - e2) == 0

	# ordering kwarg passes through to the TensorOperations parser
	f1 = zeros(Float64, s0 ⊗ s0 ← s0 ⊗ s0)
	f2 = zeros(Float64, s0 ⊗ s0 ← s0 ⊗ s0)
	@grassmann order=(2, 1) f1[1, 2; 4, 5] := a[1, 2, 3] * b[3, 4, 5]
	@tensor order=(2, 1) f2[1, 2; 4, 5] := a[1, 2, 3] * b[3, 4, 5]
	@test norm(f1 - f2) == 0
	@test norm(f1 - c1) == 0   # the contraction order does not change the result

	# a custom backend keyword is rejected (macro-expansion error, wrapped
	# in a LoadError by @eval)
	@test_throws LoadError @eval @grassmann backend = TensorOperations.BaseBackend() c[1, 2; 4, 5] := a[1, 2, 3] * b[3, 4, 5]
end

@testset "API: @grassmann vs @tensor (fermionic U-turn sign)" begin
	sf = Z2Space(0 => 1, 1 => 1)

	# a plain contraction (domain leg into codomain leg) carries no fermionic
	# sign: the two backends agree on the full Z2-graded space as well
	pa = randn(Float64, sf ⊗ sf ← sf)
	pb = randn(Float64, sf ← sf ⊗ sf)
	d1 = zeros(Float64, sf ⊗ sf ← sf ⊗ sf)
	d2 = zeros(Float64, sf ⊗ sf ← sf ⊗ sf)
	@grassmann d1[1, 2; 4, 5] := pa[1, 2, 3] * pb[3, 4, 5]
	@tensor d2[1, 2; 4, 5] := pa[1, 2, 3] * pb[3, 4, 5]
	@test norm(d1 - d2) == 0

	# a closed U-turn loop (trace) picks up a factor -1 per odd-parity loop.
	# Decompose two rank (1, 1) tensors into their parity sectors:
	#   @grassmann: even_0 * even_1 - odd_0 * odd_1   (one fermionic loop)
	#   @tensor:    even_0 * even_1 + odd_0 * odd_1   (bosonic)
	u = randn(Float64, sf ← sf)
	v = randn(Float64, sf ← sf)
	u0 = zeros(Float64, sf ← sf); copy!(block(u0, Z2Irrep(0)), block(u, Z2Irrep(0)))
	u1 = zeros(Float64, sf ← sf); copy!(block(u1, Z2Irrep(1)), block(u, Z2Irrep(1)))
	v0 = zeros(Float64, sf ← sf); copy!(block(v0, Z2Irrep(0)), block(v, Z2Irrep(0)))
	v1 = zeros(Float64, sf ← sf); copy!(block(v1, Z2Irrep(1)), block(v, Z2Irrep(1)))
	g = @grassmann u[1, 2] * v[2, 1]
	t = @tensor u[1, 2] * v[2, 1]
	te = @tensor u0[1, 2] * v0[2, 1]
	to = @tensor u1[1, 2] * v1[2, 1]
	@test g ≈ te - to atol=1.0e-12
	@test t ≈ te + to atol=1.0e-12
	@test !(g ≈ t)   # the fermionic sign is actually visible here

	# the same U-turn rule for the trace of a rank (2, 2) tensor: with all
	# odd-parity blocks zeroed there are no fermionic loops and the two
	# backends must agree exactly
	w = randn(Float64, sf ⊗ sf ← sf ⊗ sf)
	gw = @grassmann w[1, 2, 1, 2]
	tw = @tensor w[1, 2, 1, 2]
	@test !(gw ≈ tw)
	wz = deepcopy(w)
	for (f1, f2) in fusiontrees(wz)
		(isodd(f1.uncoupled[1].n) || isodd(f1.uncoupled[2].n)) && fill!(wz[f1, f2], 0)
	end
	gwz = @grassmann wz[1, 2, 1, 2]
	twz = @tensor wz[1, 2, 1, 2]
	@test gwz ≈ twz atol=0.0
end
