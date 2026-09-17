@testset "API: FockMPS & FockLattice" begin
	# Fock-space reference MPS: the vacuum and a random state
	psi = FockMPS(5)
	@test length(psi) == 5
	psi1 = FockMPS(Float64, 5)
	@test length(psi1) == 5
	r = randomfockmps(5, D=8)
	@test length(r) == 5 && norm(r) > 0

	# ordering hierarchy and the shorthand aliases
	@test M1N1() isa ImagFockOrdering && M1N1() isa FockOrdering
	@test MN() === M1N1()
	@test M1m1N1n1() isa RealFockOrdering && MmNn() === M1m1N1n1()
	@test M1N1_m1M1n1N1m2M2n2N2() isa MixedFockOrdering && MN_MmNn() === M1N1_m1M1n1N1m2M2n2N2()
	@test similargrassmannordering(M1N1) isa ImagGrassmannOrdering
	@test similargrassmannordering(M1m1N1n1()) isa RealGrassmannOrdering
	@test similargrassmannordering(MN_MmNn()) isa MixedGrassmannOrdering

	# Fock lattices and their Grassmann counterparts
	flat = FockLattice(N=3, δτ=0.1, contour=:imag)
	@test flat isa ImagFockLattice && flat isa AbstractFockLattice
	rlat_f = FockLattice(N=3, δt=0.05, contour=:real)
	@test rlat_f isa RealFockLattice
	mxfl = FockLattice(Nt=2, δt=0.05, Nτ=3, δτ=0.1, contour=:mixed)
	@test mxfl isa MixedFockLattice
	gl = similargrassmannlattice(flat)
	@test gl isa ImagGrassmannLattice
	@test gl.ordering == similargrassmannordering(flat.ordering)
	glr = similargrassmannlattice(rlat_f)
	@test glr isa RealGrassmannLattice
	@test glr.ordering == similargrassmannordering(rlat_f.ordering)
	glm = similargrassmannlattice(mxfl)
	@test glm isa MixedGrassmannLattice
	@test glm.ordering == similargrassmannordering(mxfl.ordering)
end

@testset "API: ExpNTerm & reweighting" begin
	trunc = truncdimcutoff(D=50, ϵ=1.0e-10)

	# exponential n̂ operators on the Fock lattice (no even-parity restriction)
	t = ExpNTerm(3, 1; coeff=0.5)
	@test t isa AbstractNTerm
	@test GTEMPO.positions(t) == (1, 3)   # qualified: ImpurityModelBase also exports positions
	@test t.coeff == 0.5
	@test_throws ArgumentError ExpNTerm(2, 2; coeff=1.0)

	# reweighting: multiply the phonon IF (built on the Fock lattice) into K
	rlat = GrassmannLattice(N=3, δτ=0.1, contour=:imag)
	flat = FockLattice(N=3, δτ=0.1, contour=:imag)
	pbath = bosonicbath(DiracDelta(ω=0.8, α=0.5), β=1.0)
	pcorr = correlationfunction(pbath, flat)
	mpsI_p = hybriddynamics(flat, pcorr, trunc=trunc)
	K = sysdynamics(rlat, ToulouseIM(ϵ_d=0.5), trunc=trunc)
	Kw = reweighting!(rlat, K, flat, mpsI_p, trunc=trunc)
	@test Kw isa GrassmannMPS && length(Kw) == length(rlat) && norm(Kw) > 0
	Kw2 = reweighting(rlat, sysdynamics(rlat, ToulouseIM(ϵ_d=0.5), trunc=trunc), flat, mpsI_p, trunc=trunc)
	@test distance(Kw, Kw2) / norm(Kw) < 1.0e-12
end
