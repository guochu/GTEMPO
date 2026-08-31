println("------------------------------------")
println("|          Initial state           |")
println("------------------------------------")



@testset "SIAM" begin
    for N in (0,1,10), β in (1, 10), (bands, U, μ) in ((1,0,0), (1,0,0.8), (2,1,0.5), (2,2,-1))
        lattice = GrassmannLattice(N=N, δt=0.05, contour=:real, bands=bands)
        model = AndersonIM(U, μ)

        res1 = initthermalstate(lattice, model, β)
        res2 = systhermalstate!(vacuumstate(lattice), lattice, model; β=β, trunc=notrunc())
        _normalize!(res2)
        @test distance(res1, res2) < 1e-6
    end
end

@testset "SKIM" begin
    for norb in 1:3, N in (0,1,10), β in (1, 10), (U,J,μ) in ((1,1,1), (0.7, 2.2, -0.1), (0.8, 1.1, 0.5))
        lattice = GrassmannLattice(N=N, δt=0.05, contour=:real, bands=norb*2)
        model = KanamoriIM(; U=U, J=J, μ=μ, norb=norb)

        res1 = initthermalstate(lattice, model, β)
        dis = map([100, 1000, 3000]) do n
            vac = vacuumstate(lattice)
            res2 = systhermalstate!(vac, lattice, model; β=β, trunc=notrunc(), δτ=β/n)
            distance(res1, res2)
        end
        @test all(diff(dis) .<= 1e-6)
    end
end

@testset "initfockstate2/initthermalstate2" begin
    # direct Fock -> coherent-state conversion must agree with the GTerm path
    for N in (0,1,10), β in (1, 10, Inf), (bands, U, μ) in ((1,0,0), (1,0,0.8), (2,1,0.5), (2,2,-1))
        lattice = GrassmannLattice(N=N, δt=0.05, contour=:real, bands=bands)
        model = AndersonIM(U, μ)
        @test distance(initthermalstate(lattice, model, β), initthermalstate2(lattice, model, β)) < 1e-7
    end
    for norb in 1:3, β in (1, 10, Inf), (U,J,μ) in ((1,1,1), (0.7, 2.2, -0.1), (0.8, 1.1, 0.5))
        lattice = GrassmannLattice(N=1, δt=0.05, contour=:real, bands=norb*2)
        model = KanamoriIM(; U=U, J=J, μ=μ, norb=norb)
        @test distance(initthermalstate(lattice, model, β), initthermalstate2(lattice, model, β)) < 1e-7
    end
    # generic (possibly complex) density matrices
    for bands in (1, 2, 3)
        lattice = GrassmannLattice(N=1, δt=0.05, contour=:real, bands=bands)
        d = 2^bands
        A = randn(d, d) + im*randn(d, d)
        ρ = A * A' / tr(A * A')
        @test distance(initfockstate(lattice, ρ), initfockstate2(lattice, ρ)) < 1e-7
    end
end
