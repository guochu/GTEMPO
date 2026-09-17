@testset "API: MPO" begin
	ph = z2space()
	Id = isomorphism(Float64, ph, ph)

	# a bond-2 cell encoding the identity operator
	cell = Matrix{Any}(undef, 2, 2)
	fill!(cell, 0.0)
	cell[1, 1] = 1.0
	cell[2, 2] = 1.0
	cell[1, 2] = Id
	W = SchurMPOTensor(cell)
	@test W isa SchurMPOTensor

	h = MPOHamiltonian([W, W, W])
	@test h isa MPOHamiltonian
	mpo = MPO(h)
	@test mpo isa AbstractMPO && length(mpo) == 3
	@test mpo[1] isa MPOTensor
	# the MPO constructor from plain site tensors
	@test MPO(mpo.data) isa AbstractMPO

	# the identity MPO applied to a random MPS: apply! agrees with the
	# partial-MPO multiplication, and the result is parallel to x up to the
	# fermionic sign convention of the MPO channel
	x = randomgmps(3, D=4)
	pmpo = PartialMPO([1, 2, 3], mpo.data)
	y = apply!(pmpo, deepcopy(x))
	@test y isa GrassmannMPS && length(y) == 3
	@test distance(y, pmpo * deepcopy(x)) == 0
	@test abs(dot(y, x)) / (norm(y) * norm(x)) ≈ 1.0 atol=1.0e-10

	# time-evolution MPO steppers (the cell encodes H = id, so the result
	# is e^{dt}-times the identity channel); the cell is exponentiated into
	# a sparse MPO tensor
	W1 = timeevompo(W, 0.1, WI())
	@test W1 isa GTEMPO.AbstractSparseMPOTensor
	W2 = timeevompo(W, 0.1, WII(tol=1.0e-14, maxiter=100000))
	@test W2 isa GTEMPO.AbstractSparseMPOTensor
	h2 = timeevompo(h, 0.1; alg=WII(tol=1.0e-14, maxiter=100000))
	@test h2 isa MPOHamiltonian
	u1, u2 = timeevompo(W, 0.1, ComplexStepper(WII(tol=1.0e-14, maxiter=100000)))
	@test u1 isa GTEMPO.AbstractSparseMPOTensor && u2 isa GTEMPO.AbstractSparseMPOTensor
	@test complex_stepper(0.1) == ((1 - im) * 0.1 / 2, (1 + im) * 0.1 / 2)
	@test WI() isa FirstOrderStepper && WII() isa FirstOrderStepper

	# DMRG algorithm hierarchy
	trunc = truncdimcutoff(D=32, ϵ=1.0e-12)
	@test SVDCompression(trunc) isa DMRGAlgorithm
	@test SVDCompression() isa DMRGAlgorithm
	@test SVDCompression(trunc=truncdim(16)) isa DMRGAlgorithm       # any TruncationScheme
	@test SVDCompression(trunc=truncrelerr(1.0e-10)) isa DMRGAlgorithm
	@test SVDCompression(trunc).trunc === trunc
	@test SVDCompression(truncdimcutoff(D=50, ϵ=1.0e-8)).trunc == truncdimcutoff(D=50, ϵ=1.0e-8, add_back=0)
	alg = SVDCompression(trunc)
	@test similar(alg).trunc == trunc
	@test similar(alg; trunc=truncdim(16)).trunc == truncdim(16)
	@test DMRG1(trunc=trunc) isa DMRGAlgorithm
	@test DMRG1(truncdim(16)) isa DMRGAlgorithm                      # any scheme with a bond dimension
	@test DMRG2(trunc) isa DMRGAlgorithm
	# DMRG1/DMRG2 need the bond dimension for their initial guesses:
	# schemes without D are rejected
	@test_throws MethodError DMRG1(truncerr(1.0e-10))
	@test_throws MethodError DMRG2(truncerr(1.0e-10))

	# tensor type helpers
	S = typeof(ph)
	@test randomgmps(3, D=2)[1] isa mpstensortype(S, Float64)
	@test zeros(Float64, ph ⊗ ph ← ph ⊗ ph) isa mpotensortype(S, Float64)
	@test zeros(Float64, ph ⊗ ph ← ph ⊗ ph) isa MPOTensor
	@test zeros(Float64, ph ← ph) isa MPSBondTensor
	@test bondtensortype(S, Float64) isa Type
	# boundary environments of an MPS
	y = randomgmps(3, D=4)
	@test norm(l_LL(y)) > 0 && norm(r_RR(y)) > 0
	@test norm(l_LL(y[1], y[1])) > 0
end
