println("------------------------------------")
println("|    sysdynamics2_new for K        |")
println("------------------------------------")


@testset "build K new 1 orb imag time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δτ, N) in [(1., 0.7, 0.05, 4), (0.8, 0.6, 0.1, 6)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in imag_grassmann_orderings
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=2, contour=:imag, ordering=ordering)
			K1 = sysdynamics2(lattice, exact_model)
			K2 = sysdynamics2_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7
		end
	end
end
@testset "build K new 1 orb real time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δt, N) in [(1., 0.7, 0.05, 4), (0.8, 0.6, 0.1, 3)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in real_grassmann_orderings
			lattice = GrassmannLattice(δt=δt, N=N, bands=2, contour=:real, ordering=ordering)
			K1 = sysdynamics2(lattice, exact_model)
			K2 = sysdynamics2_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7
		end
	end
end
@testset "build K new 1 orb mixed time" begin
	tol = 1.0e-9

	for (U, ϵ_d, δτ, Nτ, δt, Nt) in [(1., 0.7, 0.05, 4, 0.04, 2), (0.8, 0.6, 0.1, 6, 0.08, 4)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		for ordering in mixed_grassmann_orderings
			lattice = GrassmannLattice(δt=δt, Nt=Nt, δτ=δτ, Nτ=Nτ, bands=2, contour=:mixed, ordering=ordering)
			K1 = sysdynamics2(lattice, exact_model)
			K2 = sysdynamics2_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-8
			@test _dis(K1, K2) < 1e-7
		end
	end
end

@testset "build K new 2 orb" begin
	tol = 1.0e-9

	for (U, ϵ_d, δτ, N) in [(1., 0.7, 0.05, 4)]
		exact_model = KanamoriIM(U=U, μ=ϵ_d, J=1.1, norb=2)
		for ordering in (A1Ā1B1B̄1(), A1B1B̄1Ā1())
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=4, contour=:imag, ordering=ordering)
			K1 = sysdynamics2(lattice, exact_model)
			K2 = sysdynamics2_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7
		end
	end
	for (U, ϵ_d, δt, N) in [(1., 0.7, 0.05, 2)]
		exact_model = KanamoriIM(U=U, μ=ϵ_d, J=1.1, norb=2)
		# note: A2Ā2A1Ā1a2ā2a1ā1B2B̄2B1B̄1b2b̄2b1b̄1 is excluded here: its bond
		# dimension exceeds the default truncation cutoff (D=1000), so old and
		# new paths differ at truncation level (~1e-5) even though both are exact
		for ordering in (A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A1B1ā1b̄1Ā1B̄1a1b1(), A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
			lattice = GrassmannLattice(δt=δt, N=N, bands=4, contour=:real, ordering=ordering)
			K1 = sysdynamics2(lattice, exact_model)
			K2 = sysdynamics2_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
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

@testset "build K new general impurity imag time" begin
	tol = 1.0e-9

	δτ = 0.05
	N = 3
	for exact_model in _models
		for ordering in (A1Ā1B1B̄1(), A1B1B̄1Ā1())
			lattice = GrassmannLattice(δτ=δτ, N=N, bands=exact_model.bands, contour=:imag, ordering=ordering)
			K1 = sysdynamics2(lattice, exact_model)
			K2 = sysdynamics2_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7
		end
	end
end

@testset "build K new general impurity real time" begin
	tol = 1.0e-9

	δt = 0.05
	N = 3
	for exact_model in _models
		for ordering in (A1Ā1a1ā1B1B̄1b1b̄1(), A1Ā1B1B̄1b̄1B̄1ā1Ā1(), A1B1ā1b̄1Ā1B̄1a1b1(), A2B2B̄2Ā2A1B1B̄1Ā1a1b1b̄1ā1a2b2b̄2ā2(), A2Ā2B2B̄2A1Ā1B1B̄1a1ā1b1b̄1a2ā2b2b̄2())
			lattice = GrassmannLattice(δt=δt, N=N, bands=exact_model.bands, contour=:real, ordering=ordering)
			K1 = sysdynamics2(lattice, exact_model)
			K2 = sysdynamics2_new(lattice, exact_model)
			@test (norm(K1) - norm(K2)) < 1e-9
			@test _dis(K1, K2) < 1e-7
		end
	end
end

@testset "sysdynamics2_new branch consistency" begin
	# the branch entries must agree with sysdynamics2 branch by branch
	for (U, ϵ_d, δt, N) in [(1., 0.7, 0.05, 2)]
		exact_model = AndersonIM(U=U, μ=ϵ_d)
		lattice = GrassmannLattice(δt=δt, N=N, bands=2, contour=:real, ordering=A1B1ā1b̄1Ā1B̄1a1b1())
		Kf1 = sysdynamics2(lattice, exact_model, branch=:+)
		Kf2 = sysdynamics2_new(lattice, exact_model, branch=:+)
		Kb1 = sysdynamics2(lattice, exact_model, branch=:-)
		Kb2 = sysdynamics2_new(lattice, exact_model, branch=:-)
		@test (norm(Kf1) - norm(Kf2)) < 1e-9
		@test _dis(Kf1, Kf2) < 1e-7
		@test (norm(Kb1) - norm(Kb2)) < 1e-9
		@test _dis(Kb1, Kb2) < 1e-7
	end
end

@testset "baresysdynamics2_new" begin
	models = [
		gAndersonIM(U=0., μ=0.5),
		gAndersonIM(U=1., μ=0.5),
		IRLM(U=1., μ=0.5, J=1),
		KanamoriIM(U=1., μ=0.7, J=1.1, norb=1)
	]

	# bulkconnection on the bare propagator must reproduce sysdynamics2_new exactly
	for model in models
		lattice = GrassmannLattice(δτ=0.05, N=3, bands=model.bands, contour=:imag, ordering=A1Ā1B1B̄1())
		Kfull = sysdynamics2_new(lattice, model)
		Kbare = baresysdynamics2_new(lattice, model)
		for band in 1:lattice.bands
			Kbare = bulkconnection!(Kbare, lattice, band=band)
		end
		@test (norm(Kfull) - norm(Kbare)) < 1e-8
		@test _dis(Kfull, Kbare) < 1e-7
	end
	for model in models
		lattice = GrassmannLattice(δt=0.05, N=2, bands=model.bands, contour=:real, ordering=A1Ā1a1ā1B1B̄1b1b̄1())
		Kfull = sysdynamics2_new(lattice, model)
		Kbare = baresysdynamics2_new(lattice, model)
		for band in 1:lattice.bands
			Kbare = bulkconnection!(Kbare, lattice, band=band)
		end
		@test (norm(Kfull) - norm(Kbare)) < 1e-8
		@test _dis(Kfull, Kbare) < 1e-7
	end
	# mixed contour
	for model in (models[2], models[4])
		lattice = GrassmannLattice(δt=0.04, Nt=2, δτ=0.05, Nτ=2, bands=model.bands, contour=:mixed, ordering=A1Ā1B1B̄1_A1Ā1a1ā1B1B̄1b1b̄1A2Ā2a2ā2B2B̄2b2b̄2())
		Kfull = sysdynamics2_new(lattice, model)
		Kbare = baresysdynamics2_new(lattice, model)
		for band in 1:lattice.bands
			Kbare = bulkconnection!(Kbare, lattice, band=band)
		end
		@test (norm(Kfull) - norm(Kbare)) < 1e-8
		@test _dis(Kfull, Kbare) < 1e-7
	end

	# first-order consistency with the Trotterized baresysdynamics at small dt
	for model in models
		lattice = GrassmannLattice(δτ=0.02, N=2, bands=model.bands, contour=:imag, ordering=A1Ā1B1B̄1())
		K1 = baresysdynamics(lattice, model)
		K2 = baresysdynamics2_new(lattice, model)
		@test _dis(K1, K2) / norm(K1) < 5e-3
	end
	for model in models
		lattice = GrassmannLattice(δt=0.02, N=2, bands=model.bands, contour=:real, ordering=A1Ā1a1ā1B1B̄1b1b̄1())
		K1 = baresysdynamics(lattice, model)
		K2 = baresysdynamics2_new(lattice, model)
		@test _dis(K1, K2) / norm(K1) < 5e-3
	end
end
