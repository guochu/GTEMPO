println("------------------------------------")
println("|          Initial state           |")
println("------------------------------------")


@testset "SIAM" begin
    for N in (0,1,10), β in (1, 10), (bands, U, μ) in ((1,0,0), (1,0,0.8), (2,1,0.5), (2,2,-1))
        lattice = GrassmannLattice(N=N, δt=0.05, contour=:real, bands=bands)
        model = AndersonIM(U, μ)

        res1 = initthermalstate2(lattice, model, β)
        res2 = systhermalstate!(vacuumstate(lattice), lattice, model; β=β)
        _normalize!(res2)
        @test distance(res1, res2) < 1e-6
    end
end

@testset "SKIM" begin
    for norb in 1:3, N in (0,1,10), β in (1, 10), (U,J,μ) in ((1,1,1), (0.7, 2.2, -0.1), (0.8, 1.1, 0.5))
        lattice = GrassmannLattice(N=N, δt=0.05, contour=:real, bands=norb*2)
        model = KanamoriIM(; U=U, J=J, μ=μ, norb=norb)

        res1 = initthermalstate2(lattice, model, β)
        res2 = systhermalstate!(vacuumstate(lattice), lattice, model; β=β)
        _normalize!(res2)
        @test distance(res1, res2) < 1e-6
    end
end

@testset "initfockstate2" begin
    # generic (possibly complex) density matrices: repeated construction is
    # deterministic and two different states stay distinguishable
    for bands in (1, 2, 3)
        lattice = GrassmannLattice(N=1, δt=0.05, contour=:real, bands=bands)
        d = 2^bands
        A = randn(d, d) + im*randn(d, d)
        ρ = A * A' / tr(A * A')
        res1 = initfockstate2(lattice, FockMatrix(ρ))
        res2 = initfockstate2(lattice, FockMatrix(ρ))
        @test distance(res1, res2) < 1e-10
        B = randn(d, d) + im*randn(d, d)
        σ = B * B' / tr(B * B')
        res3 = initfockstate2(lattice, FockMatrix(σ))
        @test distance(res1, res3) > 1e-3
    end
    # large-β limit approaches the ground state projector
    for (bands, U, μ) in ((1,0,0.8), (2,1,0.5))
        lattice = GrassmannLattice(N=1, δt=0.05, contour=:real, bands=bands)
        model = AndersonIM(U, μ)
        resβ = initthermalstate2(lattice, model, 1000.)
        res0 = initthermalstate2(lattice, model, Inf)
        _normalize!(resβ)
        _normalize!(res0)
        @test distance(resβ, res0) < 1e-6
    end
end
