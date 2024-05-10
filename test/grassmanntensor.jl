println("------------------------------------")
println("|        Grassmann Tensor          |")
println("------------------------------------")


@testset "GrassmannTensor permute" begin
	s = Rep[ℤ₂](0=>5, 1=>4)
	# test 1
	m1 = TensorMap(randn, ComplexF64, s ⊗ s ⊗ s, s ⊗ s)
	m2 = permute(m1, (1,2), (3,4,5))
	for (f1, f2) in fusiontrees(m2)
		coef1 = (isodd(f2.uncoupled[1].n) && isodd(f2.uncoupled[2].n)) ? -1 : 1
		coef2 = (isodd(f2.uncoupled[1].n) && isodd(f2.uncoupled[3].n)) ? -1 : 1
        coef = coef1 * coef2 
        if coef != 1
            lmul!(coef, m2[f1, f2])
        end
	end
	m3 = permute(GrassmannTensorMap(m1), (1,2), (3,4,5))
	@test m2 == m3.data

	# test 2
	m1 = TensorMap(randn, ComplexF64, s ⊗ s, s ⊗ s ⊗ s)
	for p in [((1,2,5,4,3), ()), ((1,2,5,4), (3,)), ((1,2,5), (3,4)), ((1,), (3,4,5,2)), ((), (3,4,5,2,1))]
		m2 = permute(m1, p, copy=true)
		m3 = permute(GrassmannTensorMap(m1), p, copy=true)
		@test m2 == m3.data
	end

	# test 3
	m1 = TensorMap(randn, ComplexF64, s ⊗ s ⊗ s, s ⊗ s)
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
	m3 = permute(GrassmannTensorMap(m1), (1,3), (2,4,5))
	@test m2 == m3.data

	# test 4
	m1 = TensorMap(randn, ComplexF64, s ⊗ s, s ⊗ s)
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

	m3 = permute(GrassmannTensorMap(m1), (1,3), (2,4))
	@test m2 == m3.data
end


@testset "GrassmannTensor contract" begin
	s = Rep[ℤ₂](0=>3, 1=>4)

	# test 1
	m1 = TensorMap(randn, ComplexF64, s ⊗ s, s ⊗ s ⊗ s)
	m2 = TensorMap(randn, ComplexF64, s ⊗ s ⊗ s, s ⊗ s)

	@tensor m3[1,2;6,7] := m1[1,2,3,4,5] * m2[3,4,5,6,7]

	m1′ = GrassmannTensorMap(m1)
	m2′ = GrassmannTensorMap(m2)
	@tensor m3′[1,2;6,7] := m1′[1,2,3,4,5] * m2′[3,4,5,6,7]

	@test m3 == m3′.data

	# test 2
	m1 = TensorMap(randn, ComplexF64, s ⊗ s, s ⊗ s ⊗ s)
	m2 = TensorMap(randn, ComplexF64, s ⊗ s ⊗ s, s' ⊗ s)

	m1′ = GrassmannTensorMap(m1)
	m2′ = GrassmannTensorMap(m2)
	@tensor m3′[1,2;6,7] := m1′[1,2,3,4,5] * m2′[6,4,3,5,7]

	m22 = permute(m2′, (3,2,4), (1,5))
	@tensor m3[1,2;6,7] := m1[1,2,3,4,5] * m22.data[3,4,5,6,7]

	@test m3 == m3′.data
end
