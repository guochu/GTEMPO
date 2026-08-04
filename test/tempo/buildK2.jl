println("------------------------------------")
println("|       sysdynamics2 for K         |")
println("------------------------------------")


@testset "build K 1 orb imag time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δτ, N) in [(1., 0.7, 0.05, 4), (0.8, 0.6, 0.1, 6)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in imag_grassmann_orderings
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=2, contour=:imag, ordering=ordering)
			K1 = sysdynamics(lattice, exact_model)
			K2 = sysdynamics2(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7

			if ordering in (A1Ā1B1B̄1(), A1B1B̄1Ā1())
				K3 = sysdynamics_fast(lattice, exact_model)
				@test (norm(K2) - norm(K3)) < 1e-9
				@test _dis(K2, K3) < 1e-7
			end
		end
	end
end
@testset "build K 1 orb real time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δt, N) in [(1., 0.7, 0.05, 4), (0.8, 0.6, 0.1, 3)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in real_grassmann_orderings
			lattice = GrassmannLattice(δt=δt, N=N, bands=2, contour=:real, ordering=ordering)
			K1 = sysdynamics(lattice, exact_model)
			K2 = sysdynamics2(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7

			if ordering in (A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A1B1ā1b̄1Ā1B̄1a1b1(),  A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
				K3 = sysdynamics_fast(lattice, exact_model)
				@test (norm(K2) - norm(K3)) < 1e-9
				@test _dis(K2, K3) < 1e-7
			end
		end
	end
end
@testset "build K 1 orb mixed time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δτ, Nτ, δt, Nt) in [(1., 0.7, 0.05, 4, 0.04, 2), (0.8, 0.6, 0.1, 6, 0.08, 4)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in mixed_grassmann_orderings
			lattice = GrassmannLattice(δt=δt, Nt=Nt, δτ=δτ, Nτ=Nτ, bands=2, contour=:mixed, ordering=ordering)
			K1 = sysdynamics(lattice, exact_model)
			K2 = sysdynamics2(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7

			if ordering in (A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2(), A1Ā1B1B̄1_a1ā1A1Ā1b1b̄1B1B̄1a2ā2A2Ā2b2b̄2B2B̄2())
				K3 = sysdynamics_fast(lattice, exact_model)
				@test (norm(K2) - norm(K3)) < 1e-8
				@test _dis(K2, K3) < 1e-7
			end
		end
	end
end

@testset "build K 2 orb imag time" begin
	tol = 1.0e-9

	for (U, ϵ_d, J, norb, δτ, N) in [(1., 0.7, 0.3, 2, 0.05, 4), (0.8, 0.6, 1.1, 2, 0.1, 6)]
		exact_model = KanamoriIM(U=U, μ=ϵ_d, J=1.1, norb=norb)
		for ordering in (A1Ā1B1B̄1(), A1B1B̄1Ā1())
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=2*norb, contour=:imag, ordering=ordering)
			K1s = [accsysdynamics_fast(lattice, exact_model, scaling=s) for s in (1,3,5)]
			K2 = sysdynamics2(lattice, exact_model)
			@test issorted([_dis(K1, K2) for K1 in K1s], rev=true)

			K3 = sysdynamics_fast(lattice, exact_model)
			@test (norm(K2) - norm( K3)) < 1e-9
			@test _dis(K2, K3) < 1e-7
		end
	end
end

@testset "build K 2 orb real time" begin
	tol = 1.0e-9

	for (U, ϵ_d, J, norb, δt, N) in [(1., 0.7, 0.3, 2, 0.05, 2), (0.8, 0.6, 1.1, 2, 0.1, 3)]
		exact_model = KanamoriIM(U=U, μ=ϵ_d, J=1.1, norb=norb)
		for ordering in (A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A1B1ā1b̄1Ā1B̄1a1b1(),  A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
			lattice = GrassmannLattice(δt=δt, N=N, bands=2*norb, contour=:real, ordering=ordering)
			K1s = [accsysdynamics_fast(lattice, exact_model, scaling=s) for s in (1,3)]
			K2 = sysdynamics2(lattice, exact_model)
			@test issorted([_dis(K1, K2) for K1 in K1s], rev=true)

			K3 = sysdynamics_fast(lattice, exact_model)
			@test (norm(K2) - norm( K3)) < 1e-9
			@test _dis(K2, K3) < 1e-7
		end
	end
end
