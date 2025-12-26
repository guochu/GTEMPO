
include("util.jl")

println("------------------------------------")
println("|           DMRG tests             |")
println("------------------------------------")


@testset "MPS initializer: product state and random state" begin
	VZ2 = Z2Space(0=>1, 1=>1)
	physpaces = [VZ2 for i in 1:4]
	psi = prodmps(ComplexF64, [VZ2 for i in 1:4], [0, 1, 1, 0])
	@test length(psi) == 4
	for i in 1:length(psi)
		@test !isdual(space(psi[i], 1))
		@test !isdual(space(psi[i], 2))
		@test isdual(space(psi[i], 3))
	end
	for i in (1, length(psi)+1)
		@test !isdual(space(psi.s[i], 1))
		@test isdual(space(psi.s[i], 2))
	end
	@test psi.s[1] == one(psi.s[1])
	@test psi.s[end] == one(psi.s[end])
	@test scalartype(psi) == ComplexF64
	@test bond_dimensions(psi) == [1,1,1,1]
	@test space_l(psi) == oneunit(VZ2)
	@test space_r(psi) == oneunit(VZ2)'
	@test physical_spaces(psi) == physpaces
	@test isleftcanonical(psi)
	@test isrightcanonical(psi)
	@test !iscanonical(psi)
	psi2 = prodmps(ComplexF64, [VZ2 for i in 1:4], [1, 0, 1, 0])
	@test dot(psi, psi2) ≈ 0.

	VSU2 = Rep[U₁×SU₂]((-0.5, 0)=>1, (0.5, 0)=>1, (0, 0.5)=>1)
	physpaces = [VSU2 for i in 1:4]
	psi = prodmps(Float64, physpaces, [(0, 0.5), (0.5, 0), (0, 0.5), (-0.5, 0)])
	for i in 1:length(psi)
		@test !isdual(space(psi[i], 1))
		@test !isdual(space(psi[i], 2))
		@test isdual(space(psi[i], 3))
	end
	for i in (1, length(psi)+1)
		@test !isdual(space(psi.s[i], 1))
		@test isdual(space(psi.s[i], 2))
	end
	@test psi.s[1] == one(psi.s[1])
	@test psi.s[end] == one(psi.s[end])
	@test scalartype(psi) == Float64
	@test bond_dimension(psi) == 2
	@test space_l(psi) == oneunit(VSU2)
	@test space_r(psi) == oneunit(VSU2)'
	@test physical_spaces(psi) == physpaces

	# nontrivial sector
	right = Rep[U₁×SU₂]((0, 0.5)=>1)
	psi = prodmps(Float64, [VSU2 for i in 1:3], [(0.5, 0), (-0.5, 0.), (0, 0.5)], right=right)
	@test space_r(psi) == right'
	@test !isleftcanonical(psi)
	@test !isrightcanonical(psi)
	@test !iscanonical(psi)
end

@testset "MPO initializer: product operator" begin
	# u1 symmetry
	p = spin_site_ops_u1()
	sp, sm, z = p["+"], p["-"], p["z"]
	ph = space(z, 1)
	physpaces = [ph for i in 1:4]
	h1 = prodmpo(Float64, physpaces, [2], [sp])
	@test length(h1) == 4
	@test scalartype(h1) == Float64
	for i in 1:length(h1)
		@test !isdual(space(h1[i], 1))
		@test !isdual(space(h1[i], 2))
		@test isdual(space(h1[i], 3))
		@test isdual(space(h1[i], 4))
	end	
	@test space_l(h1) == oneunit(space_l(h1))
	@test space_r(h1) == space_r(sp)

	h1 = prodmpo(Float64, physpaces, [1, 3], [sp, sp'])
	@test space_l(h1) == oneunit(space_l(h1))
	@test space_r(h1)' == oneunit(space_r(h1))

end

@testset "Exponential expansion    " begin
	L = 100
	atol = 1.0e-5
	for alpha in (-2, -2.5, -3)
		xdata = [convert(Float64, i) for i in 1:L]
		ydata = [1.3 * x^alpha for x in xdata]
		xs1, lambdas1 = exponential_expansion(ydata, PronyExpansion(n=20,tol=atol))
		@test expansion_error(ydata, xs1, lambdas1) < atol
		xs2, lambdas2 = exponential_expansion(ydata, DeterminedPronyExpansion(n=20,tol=atol))
		@test expansion_error(ydata, xs2, lambdas2) < atol
		xs3, lambdas3 = exponential_expansion(ydata, PronyExpansion2(n=20,atol=atol))
		@test expansion_error(ydata, xs3, lambdas3) < atol
		# xs2, lambdas2 = exponential_expansion(ydata, LsqExpansion(atol=atol))
		# @test expansion_error(ydata, xs2, lambdas2) < atol
	end
	L = 500
	xdata = [convert(Float64, i) for i in 1:L]
	ydata = [1.3 * 0.7^x + 0.7 * 0.5^x - 1.1 * 0.8^x + 0.1*0.95^x for x in xdata]
	xs1, lambdas1 = exponential_expansion(ydata, PronyExpansion(n=20,tol=atol))
	@test expansion_error(ydata, xs1, lambdas1) < atol
	for stepsize in 2:6
		xs2, lambdas2 = exponential_expansion(ydata[1:stepsize:L], PronyExpansion(n=20, stepsize=stepsize, tol=atol))
		@test expansion_error(ydata, xs2, lambdas2) < atol
	end
	L = 500
	xdata = [convert(Float64, i) for i in 1:L]
	ydata = [(1.3+0.2im) * (0.7+0.3im)^x + (0.7+1.1im) * (0.5+0.1im)^x - (1.1-0.3im) * (0.8-0.5im)^x + (0.1-0.2)*(0.95+0.2im)^x for x in xdata]
	for alg in (PronyExpansion(n=20,tol=atol), PronyExpansion2(n=20,atol=atol), LsqExpansion2(n=20,atol=atol))
		xs, lambdas = exponential_expansion(ydata, alg)
		@test expansion_error(ydata, xs, lambdas) < atol
	end
end
