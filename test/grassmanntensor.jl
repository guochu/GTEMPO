println("------------------------------------")
println("|        Grassmann Tensor          |")
println("------------------------------------")


@testset "GrassmannTensor permute" begin
	s = Z2Space(0=>5, 1=>4)
	# test 1
	m1 = randn(ComplexF64, s ⊗ s ⊗ s, s ⊗ s)
	m2 = permute(m1, (1,2), (3,4,5))
	for (f1, f2) in fusiontrees(m2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
        coef = coef1 * coef2
        if coef != 1
            lmul!(coef, m2[f1, f2])
        end
	end
	@grassmann m3[1,2; 3,4,5] := m1[1,2,3,4,5]
	@test m2 == m3

	# test 2
	m1 = randn(ComplexF64, s ⊗ s, s ⊗ s ⊗ s)
	# these permutations carry no fermionic sign; the fermionic permute is
	# expressed through @grassmann (the empty-codomain permutation of the same
	# sign class cannot be written in macro syntax)
	m2 = permute(m1, ((1,2,5,4,3), ()), copy=true)
	@grassmann m3[1,2,5,4,3;] := m1[1,2,3,4,5]
	@test m2 == m3

	m2 = permute(m1, ((1,2,5,4), (3,)), copy=true)
	@grassmann m3[1,2,5,4; 3] := m1[1,2,3,4,5]
	@test m2 == m3

	m2 = permute(m1, ((1,2,5), (3,4)), copy=true)
	@grassmann m3[1,2,5; 3,4] := m1[1,2,3,4,5]
	@test m2 == m3

	m2 = permute(m1, ((1,), (3,4,5,2)), copy=true)
	@grassmann m3[1; 3 4 5 2] := m1[1,2,3,4,5]
	@test m2 == m3

	# test 3
	m1 = randn(ComplexF64, s ⊗ s ⊗ s, s ⊗ s)
	m2 = permute(m1, (1,3), (2,4,5))
	for (f1, f2) in fusiontrees(m2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
		coef3 = (isodd(f2.uncoupled[1].n) && isodd(f1.uncoupled[2].n)) ? -1 : 1

        coef = coef1 * coef2 * coef3
        if coef != 1
            lmul!(coef, m2[f1, f2])
        end
	end
	@grassmann m3[1,3; 2,4,5] := m1[1,2,3,4,5]
	@test m2 == m3

	# test 4
	m1 = randn(ComplexF64, s ⊗ s, s ⊗ s)
	m2 = permute(m1, (1,3), (2,4))
	for (f1, f2) in fusiontrees(m2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[1].n)) ? -1 : 1
		coef3 = (isodd(f1.uncoupled[2].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1

        coef = coef1 * coef2 * coef3
        if coef != 1
            lmul!(coef, m2[f1, f2])
        end
	end

	@grassmann m3[1,3; 2,4] := m1[1,2,3,4]
	@test m2 == m3
end


@testset "GrassmannTensor contract" begin
	s = Z2Space(0=>3, 1=>4)

	# test 1
	m1 = randn(ComplexF64, s ⊗ s, s ⊗ s ⊗ s)
	m2 = randn(ComplexF64, s ⊗ s ⊗ s, s ⊗ s)

	@tensor m3[1,2;6,7] := m1[1,2,3,4,5] * m2[3,4,5,6,7]

	@grassmann m3′[1,2;6,7] := m1[1,2,3,4,5] * m2[3,4,5,6,7]

	@test m3 == m3′

	# test 2
	m1 = randn(ComplexF64, s ⊗ s, s ⊗ s ⊗ s)
	m2 = randn(ComplexF64, s ⊗ s ⊗ s, s' ⊗ s)

	@grassmann m3′[1,2;6,7] := m1[1,2,3,4,5] * m2[6,4,3,5,7]

	@grassmann m22[3,2,4; 1,5] := m2[1,2,3,4,5]
	@tensor m3[1,2;6,7] := m1[1,2,3,4,5] * m22[3,4,5,6,7]

	@test m3 == m3′
end

@testset "@grassmann macro" begin
	s = Z2Space(0=>3, 1=>4)
	m1 = randn(ComplexF64, s ⊗ s, s ⊗ s ⊗ s)
	m2 = randn(ComplexF64, s ⊗ s ⊗ s, s ⊗ s)

	# definition of a new tensor (this contraction pattern carries no fermionic sign)
	@grassmann c[1,2;6,7] := m1[1,2,3,4,5] * m2[3,4,5,6,7]
	@tensor c′[1,2;6,7] := m1[1,2,3,4,5] * m2[3,4,5,6,7]
	@test c == c′
	@test c isa TensorMap

	# even-only tensors: fermionic and bosonic contractions coincide,
	# so all macro machinery (assignment forms, scalars, ...) can be
	# checked against the plain @tensor results
	se = Z2Space(0=>3)
	e1 = randn(ComplexF64, se ⊗ se, se ⊗ se ⊗ se)
	e2 = randn(ComplexF64, se ⊗ se ⊗ se, se ⊗ se)

	# in-place assignment forms
	c2 = randn(ComplexF64, se ⊗ se, se ⊗ se)
	r2 = copy(c2)
	@grassmann c2[1,2;6,7] = e1[1,2,3,4,5] * e2[3,4,5,6,7]
	@tensor r2[1,2;6,7] = e1[1,2,3,4,5] * e2[3,4,5,6,7]
	@test c2 == r2

	c3 = randn(ComplexF64, se ⊗ se, se ⊗ se)
	r3 = copy(c3)
	@grassmann c3[1,2;6,7] += e1[1,2,3,4,5] * e2[3,4,5,6,7]
	@tensor r3[1,2;6,7] += e1[1,2,3,4,5] * e2[3,4,5,6,7]
	@test c3 == r3

	# scalar output
	@grassmann s1 = e1[1,2,3,4,5] * e2[3,4,5,1,2]
	@tensor s2 = e1[1,2,3,4,5] * e2[3,4,5,1,2]
	@test s1 ≈ s2

	# conjugation (svdmult-style pattern: the contracted indices of the conj'd
	# tensor are its domain indices)
	a4 = randn(ComplexF64, se ⊗ se, se ⊗ se)
	l = randn(ComplexF64, se, se ⊗ se)
	@grassmann w[1,2;5] := a4[1,2,3,4] * conj(l[5,3,4])
	@tensor w′[1,2;5] := a4[1,2,3,4] * conj(l[5,3,4])
	@test w == w′

	# tensor object given by an arbitrary expression
	ms = [m1]
	@grassmann c4[1,2;6,7] := ms[1][1,2,3,4,5] * m2[3,4,5,6,7]
	@test c4 == c′

	# inside a comprehension
	cs = [@grassmann tmp[1,2;6,7] := m1[1,2,3,4,5] * m2[3,4,5,6,7] for _ in 1:2]
	@test cs[1] == c′
	@test cs[2] == c′

	# keyword argument form (same keywords as @tensor)
	@grassmann contractcheck=true c6[1,2;6,7] := e1[1,2,3,4,5] * e2[3,4,5,6,7]
	@tensor c6′[1,2;6,7] := e1[1,2,3,4,5] * e2[3,4,5,6,7]
	@test c6 == c6′
end

# regression test for the removal of the g_rightorth / g_stable_tsvd wrappers:
# on a single site tensor the bracketing fermionic permutes of a right sweep
# cancel pairwise, so factorizing the bosonically permuted tensor and restoring
# it bosonically stores exactly the same tensors as the old fermionic path
@testset "single-site factorizations: bosonic == fermionic" begin
	Dl = Z2Space(0=>3, 1=>2)
	ph = Z2Space(0=>4, 1=>3)
	Dr = Z2Space(0=>2, 1=>3)
	for _ in 1:3
		t = randn(ComplexF64, Dl ⊗ ph, Dr)

		# the two fermionic permutes cancel: they are involutive here
		@grassmann tf[1; 2 3] := t[1,2,3]
		@grassmann t2[1 2; 3] := tf[1,2,3]
		@test t2 ≈ t

		# LQ split: bosonic path stores the same tensors as the fermionic path
		@grassmann tf1[1; 2 3] := t[1,2,3]
		lf, qf = TK.rightorth!(tf1; alg=LQ())
		l, q = TK.rightorth!(permute(t, (1,), (2, 3); copy=true); alg=LQ())
		@test l ≈ lf
		@grassmann qf2[1 2; 3] := qf[1,2,3]
		@test permute(q, (1, 2), (3,); copy=true) ≈ qf2
		lqf = lf * qf
		@grassmann lqf2[1 2; 3] := lqf[1,2,3]
		@test lqf2 ≈ t

		# SVD split
		@grassmann tf2[1; 2 3] := t[1,2,3]
		uf, sf, vf, _ = tsvd(tf2; alg=SDD())
		u, s, v, _ = tsvd(permute(t, (1,), (2, 3); copy=true); alg=SDD())
		@test s ≈ sf
		@grassmann vf2[1 2; 3] := vf[1,2,3]
		@test permute(v, (1, 2), (3,); copy=true) ≈ vf2
		usv = uf * sf * vf
		@grassmann usv2[1 2; 3] := usv[1,2,3]
		@test usv2 ≈ t
	end
end
