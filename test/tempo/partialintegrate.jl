println("------------------------------------")
println("|        Partial Integrate         |")
println("------------------------------------")


function _normalize!(psi::GrassmannMPS)
    alg = Orthogonalize(SVD(), normalize=true)
    open("/dev/null", "w") do devnull
        redirect_stderr(devnull) do
            canonicalize!(psi, alg=alg)
        end
    end

end
function _dis(mps1, mps2)
    mps1, mps2 = deepcopy(mps1), deepcopy(mps2)
    _normalize!(mps1)
    _normalize!(mps2)
    return distance(mps1, mps2)
end

@testset "multiple GrassmannMPS multiplication" begin
	L = 8
    N = 6
    δ = 1e-7
	for T in (Float64, ComplexF64)
        psi = randomgmps(T, L, D=4)
        res0 = psi
        xs = [psi]
        for i in 2:N
            psi = randomgmps(T, L, D=4)
            push!(xs, psi)
            res0 = mult(res0, psi)
            res1 = my_mult(xs...)
            res2 = my_mult2(xs...)
            @test distance(res0, res1) < δ
            @test distance(res0, res2) < δ
        end
    end
end

@testset "multiple GrassmannMPS partialintegrate" begin
    lattice = GrassmannLattice(N=1, δτ=0.1, bands=3, contour=:imag, ordering=A1Ā1B1B̄1())
    lattice2 = similar(lattice, bands=2)
    L = length(lattice)
	trunc1 = truncdimcutoff(D=1024, ϵ=1.0e-10)
	trunc2 = truncdimcutoff(D=500, ϵ=1.0e-10)
    alg1 = SVDCompression(trunc2)
    alg2 = DMRGMult1(trunc2, initguess=:svd, maxiter=5, verbosity=0)
    algs = [alg1, alg2]

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

