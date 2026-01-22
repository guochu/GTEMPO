println("------------------------------------")
println("|    Partial Integrate in CUDA     |")
println("------------------------------------")



function _dis(mps1, mps2)
    mps1, mps2 = deepcopy(mps1), deepcopy(mps2)
    _normalize!(mps1)
    _normalize!(mps2)
    return distance(mps1, mps2)
end


@testset "multiple GrassmannMPS partialintegrate" begin
    lattice = GrassmannLattice(N=1, δτ=0.1, bands=3, contour=:imag, ordering=A1Ā1B1B̄1())
    lattice2 = similar(lattice, bands=2)
    L = length(lattice)
	trunc1 = truncdimcutoff(D=1024, ϵ=1.0e-10)
	trunc2 = truncdimcutoff(D=500, ϵ=1.0e-10)
    alg1 = CuSVDCompression(trunc2)
    alg2 = CuDMRGMult1(trunc2, initguess=:svd, maxiter=5, verbosity=0)
    algs = [alg1, alg2]
    algs = [alg2,]

    branchs = (:τ,)
    bands = [(1,), (2,), (3,), (1,2), (1,3), (2,3)]

    δ = 1e-7
	for T in (Float64, ComplexF64)
        N = 5
        psi = randomgmps(T, L, D=4)
        _normalize!(psi)
        setscaling!(psi, 0.98)
        res0 = psi
        xs = [psi]
    
        for i in 2:N
            psi = randomgmps(T, L, D=4)
            _normalize!(psi)
            setscaling!(psi, 0.98)
            push!(xs, psi)
            res0 = mult(res0, psi, trunc=trunc1)

            for j in eachindex(bands)
                res1 = integratebands(lattice, res0, bands[j])
                res2 = partialintegrate(lattice, alg1, xs...; branchs=branchs, bands=bands[j])
                res3 = partialintegrate(lattice, alg2, xs...; branchs=branchs, bands=bands[j])
                @test distance(res1, res2) < δ
                @test _dis(res1, res3) < δ    
            end
        end
    end
end

