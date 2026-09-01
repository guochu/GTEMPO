println("------------------------------------")
println("|          Initial state           |")
println("------------------------------------")


@testset "SIAM" begin
	# the analytical GTerm solution for AndersonIM must coincide with the
	# generic Fock-space construction on an equivalent ImpurityHamiltonian
	# (same operator, including the normalization)
	for N in (0,1,10), β in (1, 10), (bands, U, μ) in ((1,0,0), (1,0,0.8), (2,1,0.5), (2,2,-1))
		lattice = GrassmannLattice(N=N, δt=0.05, contour=:real, bands=bands)
		model = AndersonIM(U, μ)

		res1 = systhermalstate!(vacuumstate(lattice), lattice, model; β=β)
		res2 = systhermalstate!(vacuumstate(lattice), lattice, gAndersonIM(U=U, μ=μ); β=β)
		@test _dis(res1, res2) < 1e-6
		@test abs(norm(res1) - norm(res2)) / norm(res2) < 1e-6
	end
	# β = Inf: ground state projector, analytical dispatch included
	for (bands, U, μ) in ((1,0,0.8), (2,1,0.5), (2,1,-0.5))
		lattice = GrassmannLattice(N=1, δt=0.05, contour=:real, bands=bands)
		model = AndersonIM(U, μ)
		res1 = systhermalstate!(vacuumstate(lattice), lattice, model; β=Inf)
		res2 = systhermalstate!(vacuumstate(lattice), lattice, gAndersonIM(U=U, μ=μ); β=Inf)
		@test _dis(res1, res2) < 1e-6
	end
end

@testset "SKIM" begin
	for norb in 1:3, N in (0,1,10), β in (1, 10), (U,J,μ) in ((1,1,1), (0.7, 2.2, -0.1), (0.8, 1.1, 0.5))
		lattice = GrassmannLattice(N=N, δt=0.05, contour=:real, bands=norb*2)
		model = KanamoriIM(; U=U, J=J, μ=μ, norb=norb)

		res1 = systhermalstate(lattice, model; β=β)
		res2 = systhermalstate!(vacuumstate(lattice), lattice, model; β=β)
		@test _dis(res1, res2) < 1e-6
	end
end

@testset "sysinitialstate" begin
	# generic (possibly complex) density matrices: repeated construction is
	# deterministic and two different states stay distinguishable
	for bands in (1, 2, 3)
		lattice = GrassmannLattice(N=1, δt=0.05, contour=:real, bands=bands)
		d = 2^bands
		A = randn(d, d) + im*randn(d, d)
		ρ = A * A' / tr(A * A')
		res1 = sysinitialstate(lattice, FockMatrix(ρ))
		res2 = sysinitialstate!(vacuumstate(lattice), lattice, FockMatrix(ρ))
		@test distance(res1, res2) < 1e-10
		B = randn(d, d) + im*randn(d, d)
		σ = B * B' / tr(B * B')
		res3 = sysinitialstate(lattice, FockMatrix(σ))
		@test distance(res1, res3) > 1e-3
	end
	# large-β limit approaches the ground state projector
	for (bands, U, μ) in ((1,0,0.8), (2,1,0.5))
		lattice = GrassmannLattice(N=1, δt=0.05, contour=:real, bands=bands)
		model = AndersonIM(U, μ)
		resβ = systhermalstate(lattice, model; β=1000.)
		res0 = systhermalstate(lattice, model; β=Inf)
		_normalize!(resβ)
		_normalize!(res0)
		@test distance(resβ, res0) < 1e-6
	end
end
