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
	m3 = g_permute(m1, (1,2), (3,4,5))
	@test m2 == m3

	# test 2
	m1 = randn(ComplexF64, s ⊗ s, s ⊗ s ⊗ s)
	for p in [((1,2,5,4,3), ()), ((1,2,5,4), (3,)), ((1,2,5), (3,4)), ((1,), (3,4,5,2)), ((), (3,4,5,2,1))]
		m2 = permute(m1, p, copy=true)
		m3 = g_permute(m1, p, copy=true)
		@test m2 == m3
	end

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
	m3 = g_permute(m1, (1,3), (2,4,5))
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

	m3 = g_permute(m1, (1,3), (2,4))
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

	m22 = g_permute(m2, (3,2,4), (1,5))
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
