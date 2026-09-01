println("------------------------------------")
println("|       sysdynamics2 for K         |")
println("------------------------------------")


@testset "build K 1 orb imag time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δτ, N) in [(1., 0.7, 0.05, 4), (0.8, 0.6, 0.1, 6)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in imag_grassmann_orderings
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=2, contour=:imag, ordering=ordering)
			K1 = sysdynamics2_new(lattice, exact_model)
			K2 = sysdynamics2_fast_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7
		end
	end
end
@testset "build K 1 orb real time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δt, N) in [(1., 0.7, 0.05, 4), (0.8, 0.6, 0.1, 3)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in real_grassmann_orderings
			lattice = GrassmannLattice(δt=δt, N=N, bands=2, contour=:real, ordering=ordering)
			K1 = sysdynamics2_new(lattice, exact_model)
			K2 = sysdynamics2_fast_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7
		end
	end
end
@testset "build K 1 orb mixed time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δτ, Nτ, δt, Nt) in [(1., 0.7, 0.05, 4, 0.04, 2), (0.8, 0.6, 0.1, 6, 0.08, 4)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in mixed_grassmann_orderings
			lattice = GrassmannLattice(δt=δt, Nt=Nt, δτ=δτ, Nτ=Nτ, bands=2, contour=:mixed, ordering=ordering)
			K1 = sysdynamics2_new(lattice, exact_model)
			K2 = sysdynamics2_fast_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-8
			@test _dis(K1, K2) < 1e-7
		end
	end
end

@testset "build K 2 orb imag time" begin
	tol = 1.0e-9

	for (U, ϵ_d, J, norb, δτ, N) in [(1., 0.7, 0.3, 2, 0.05, 4), (0.8, 0.6, 1.1, 2, 0.1, 6)]
		exact_model = KanamoriIM(U=U, μ=ϵ_d, J=1.1, norb=norb)
		for ordering in (A1Ā1B1B̄1(), A1B1B̄1Ā1())
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=2*norb, contour=:imag, ordering=ordering)
			K1 = sysdynamics2_new(lattice, exact_model)
			K2 = sysdynamics2_fast_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-8
			@test _dis(K1, K2) < 1e-7
		end
	end
end

@testset "build K 2 orb real time" begin
	tol = 1.0e-9

	for (U, ϵ_d, J, norb, δt, N) in [(1., 0.7, 0.3, 2, 0.05, 2), (0.8, 0.6, 1.1, 2, 0.1, 3)]
		exact_model = KanamoriIM(U=U, μ=ϵ_d, J=1.1, norb=norb)
		for ordering in (A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A1B1ā1b̄1Ā1B̄1a1b1(),  A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
			lattice = GrassmannLattice(δt=δt, N=N, bands=2*norb, contour=:real, ordering=ordering)
			K1 = sysdynamics2_new(lattice, exact_model)
			K2 = sysdynamics2_fast_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-8
			@test _dis(K1, K2) < 1e-7
		end
	end
end



_models = [
	gAndersonIM(U=0., μ=0.5),
	gAndersonIM(U=1., μ=0.5),
	IRLM(U=0., μ=0.5, J=1),
	IRLM(U=1., μ=0.5, J=1),
	KanamoriIM(U=1., μ=0.7, J=1.1, norb=1),
	KanamoriIM(U=1., μ=0.7, J=1.1, norb=2)
]

@testset "build K general impurity imag time" begin
	tol = 1.0e-9

	δτ = 0.05
	N = 3
	for exact_model in _models
		for ordering in (A1Ā1B1B̄1(), A1B1B̄1Ā1())
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=exact_model.bands, contour=:imag, ordering=ordering)
			K1 = sysdynamics2_new(lattice, exact_model)
			K2 = sysdynamics2_fast_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-8
			@test _dis(K1, K2) < 1e-7
		end
	end
end

@testset "build K general impurity real time" begin
	tol = 1.0e-9

	δt = 0.05
	N = 3
	for exact_model in _models
		for ordering in (A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A1B1ā1b̄1Ā1B̄1a1b1(),  A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
			lattice = GrassmannLattice(δt=δt, N=N, bands=exact_model.bands, contour=:real, ordering=ordering)
			K1 = sysdynamics2_new(lattice, exact_model)
			K2 = sysdynamics2_fast_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-8
			@test _dis(K1, K2) < 1e-7
		end
	end
end
