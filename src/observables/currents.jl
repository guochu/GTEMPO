"""
    electriccurrent(lattice::RealGrassmannLattice1Order, corr, k::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                         alg::IntegrationAlgorithm=ExactIntegrate(),
                         Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1, max_range::Int=10000)

Return the electric current at time step k, which is calculated as the sum
of a series of Green's functions
"""
function electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                         alg::IntegrationAlgorithm=ExactIntegrate(),
                         Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    for j in 1:k-1
        curr += η⁺⁺[k, j] * Gt(lattice, k, j, A, B...; b1=:+, b2=:+, Z=Z, band=band, alg=alg, c1=true, c2=false)
        curr += η⁺⁻[k, j] * Gt(lattice, k, j, A, B...; b1=:+, b2=:-, Z=Z, band=band, alg=alg, c1=true, c2=false)
    end
    return 2 * curr / lattice.δt
end
electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; band::Int=1,
                alg::IntegrationAlgorithm=ExactIntegrate(), 
                Z::Number = integrate(lattice, A, B..., alg=alg)) = [electriccurrent(lattice, corr, k, A, B...; alg=alg, Z=Z, band=band) for k in 2:lattice.k]

"""
    electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr, k::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                         alg::IntegrationAlgorithm=ExactIntegrate(),
                         Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1, max_range::Int=10000)

The same as electriccurrent, but use a very efficient algorithm which directly 
builds the current operator as an MPO of bond dimension 2
"""
function electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::GrassmannMPS, B::Vararg{GrassmannMPS}; band::Int=1, 
                            alg::IntegrationAlgorithm=ExactIntegrate(), 
                            Z::Number = integrate(lattice, A, B..., alg=alg))
    mpo = build_current_mpo(lattice, corr, k, band)
    A2 = mpo * A
    curr = integrate(lattice, A2, B..., alg=alg) / Z
    return 2 * curr / lattice.δt
end
electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, A::GrassmannMPS, B::Vararg{GrassmannMPS}; band::Int=1,
                   alg::IntegrationAlgorithm=ExactIntegrate(),
                   Z::Number = integrate(lattice, A, B..., alg=alg)) = [electriccurrent_fast(lattice, corr, k, A, B...; alg=alg, Z=Z, band=band) for k in 2:lattice.k]



function heatcorrelationfunction(bath::AbstractFermionicBath, lattice::RealGrassmannLattice)
    bath2 = similar(bath, _mult_w(bath.spectrum))
    corr = correlationfunction(bath2, lattice)
    return corr
end

_mult_w(x::AbstractSpectrumFunction) = similar(x, w -> w * x.f(w)) 
_mult_w(x::DiracDelta) = similar(x, α=x.ω*x.α)


# """
#     heatcurrent_fast(lattice::RealGrassmannLattice1Order, bath::AbstractFermionicBath, args...; kwargs...)

# The calculation of heat current is very similar to electric current, 
# the only change one needs to make is to replace the spectrum function
# as J(w) to w -> w * J(w)
# """
heatcurrent_fast(lattice::RealGrassmannLattice1Order, bath::AbstractFermionicBath, args...; kwargs...) = electriccurrent_fast(
                    lattice, heatcorrelationfunction(bath, lattice), args...; kwargs...)




function electriccurrent(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::GrassmannMPS, B::Vararg{GrassmannMPS}; 
                            alg::IntegrationAlgorithm=ExactIntegrate(), 
                            Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    @assert  lattice.N <= div(size(η⁺⁺,1)-1, 2) 
    k = lattice.k
    curr -= η⁺⁺[2*k-1, 1] * Gt(lattice, k, 1, A, B...; b1=:-, b2=:+, Z=Z, band=band, alg=alg, c1=true, c2=false)
    curr -= η⁺⁻[2*k-1, 1] * Gt(lattice, k, 1, A, B..., b1=:-, b2=:-, Z=Z, band=band, alg=alg, c1=true, c2=false)
    for j in 2:k-1
        curr -= (η⁺⁺[2*k-1, 2*j-2] + η⁺⁺[2*k-1, 2*j-1]) * Gt(lattice, k, j, A, B..., b1=:-, b2=:+, Z=Z, band=band, alg=alg, c1=true, c2=false)
        curr -= (η⁺⁻[2*k-1, 2*j-2] + η⁺⁻[2*k-1, 2*j-1])  * Gt(lattice, k, j, A, B..., b1=:-, b2=:-, Z=Z, band=band, alg=alg, c1=true, c2=false)
    end
    curr -= η⁺⁺[2*k-1, 2*k-2] * Gt(lattice, k, k, A, B..., b1=:-, b2=:+, Z=Z, band=band, alg=alg, c1=true, c2=false)
    curr -= η⁺⁻[2*k-1, 2*k-2] * Gt(lattice, k, k, A, B..., b1=:-, b2=:-, Z=Z, band=band, alg=alg, c1=true, c2=false)
    return 2 * curr / (0.5*lattice.δt)
end

function electriccurrent_fast(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::GrassmannMPS, B::Vararg{GrassmannMPS}; 
                            alg::IntegrationAlgorithm=ExactIntegrate(), 
                            Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1)
    k = lattice.k
    mpo = build_current_mpo(lattice, corr, k, band)
    A2 = mpo * A
    curr = integrate(lattice, A2, B..., alg=alg) / Z

    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    curr -= η⁺⁺[2*k-1, 2*k-2] * Gt(lattice, k, k, A, B..., b1=:-, b2=:+, Z=Z, band=band, alg=alg, c1=true, c2=false)
    curr -= η⁺⁻[2*k-1, 2*k-2] * Gt(lattice, k, k, A, B..., b1=:-, b2=:-, Z=Z, band=band, alg=alg, c1=true, c2=false)
    return 2 * curr / (0.5*lattice.δt)
end

# building the current operator as MPO
build_current_mpo(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, band::Int) = build_current_mpo_1order(lattice, corr, k, band)

function build_current_mpo_1order(lattice::RealGrassmannLattice, corr::RealCorrelationFunction, k::Int, band::Int)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    row = index(lattice, k, conj=true, branch=:+, band=band)
    cols = Int[]
    coefs = ComplexF64[]
    for j in k-1:-1:1
        col1 = index(lattice, j, conj=false, branch=:+, band=band)
        push!(cols, col1)
        push!(coefs, η⁺⁺[k, j])
        col2 = index(lattice, j, conj=false, branch=:-, band=band)
        push!(cols, col2)
        push!(coefs, η⁺⁻[k, j])
    end
    return two_body_mpo_row(row, cols, coefs)
end

function build_current_mpo(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, k::Int, band::Int)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    row = index(lattice, k, conj=true, branch=:-, band=band)
    cols = Int[]
    coefs = ComplexF64[]
    # the last two terms are not includes
    for j in k-1:-1:2
        col1 = index(lattice, j, conj=false, branch=:+, band=band)
        push!(cols, col1)
        push!(coefs, -(η⁺⁺[2*k-1, 2*j-2] + η⁺⁺[2*k-1, 2*j-1]) )
        col2 = index(lattice, j, conj=false, branch=:-, band=band)
        push!(cols, col2)
        push!(coefs, -((η⁺⁻[2*k-1, 2*j-2] + η⁺⁻[2*k-1, 2*j-1])))
    end
    j = 1
    col1 = index(lattice, j, conj=false, branch=:+, band=band)
    push!(cols, col1)
    push!(coefs, -η⁺⁺[2*k-1, 2*j-1]) 
    col2 = index(lattice, j, conj=false, branch=:-, band=band)
    push!(cols, col2)
    push!(coefs, -η⁺⁻[2*k-1, 2*j-1])
   
    return two_body_mpo_row(row, cols, coefs)
end

function two_body_mpo_row(row::Int, cols::Vector{Int}, coefs::Vector{<:Number}) 
    @assert length(cols) == length(coefs) >= 1
    p = sortperm(cols)
    cols = cols[p]
    coefs = coefs[p]
    # @assert all(x -> row < x, cols) 
    I2 = one(JW)

    virtual = isomorphism(eltype(coefs), Rep[ℤ₂](1=>1), Rep[ℤ₂](1=>1))
    pspace = grassmannpspace()
    T = scalartype(virtual)
    @tensor m22JW[1,3;2,4] := virtual[1,2] * JW[3,4] 
    if row < cols[1]
        tmp = Matrix{Any}(undef, 1, 2)
        tmp[1, 1] = zero(I2)
        tmp[1, 2] = one(eltype(coefs)) * σ₊
        mpoj = SparseMPOTensor(tmp, T, pspace)
        data = [mpoj]
        for i in row+1 : cols[end]-1
            tmp = Matrix{Any}(undef, 2, 2)
            tmp[1,1] = I2
            tmp[2,2] = m22JW
            tmp[1,2] = 0.
            pos = findfirst(x->x==i, cols)
            if isnothing(pos)
                tmp[2,1] = 0.
            else
                tmp[2,1] = adjoint(σ₋) * coefs[pos]
            end
            push!(data, SparseMPOTensor(tmp, T, pspace))
        end
        tmp = Matrix{Any}(undef, 2, 1)
        tmp[1,1] = I2
        tmp[2,1] = adjoint(σ₋) * coefs[end]
        push!(data, SparseMPOTensor(tmp, T, pspace))
        positions = collect(row:cols[end])
    elseif row > cols[end]
        tmp = Matrix{Any}(undef, 1, 2)
        tmp[1, 1] = I2
        tmp[1, 2] = (-coefs[1]) * σ₊
        mpoj = SparseMPOTensor(tmp, T, pspace)
        data = [mpoj]
        for i in cols[1]+1:row-1
            tmp = Matrix{Any}(undef, 2, 2)
            tmp[1,1] = I2
            tmp[2,2] = m22JW
            tmp[2,1] = 0.
            pos = findfirst(x->x==i, cols)
            if isnothing(pos)
                tmp[1,2] = 0.
            else
                tmp[1,2] = (-coefs[pos]) * σ₊
            end
            push!(data, SparseMPOTensor(tmp, T, pspace))
        end
        tmp = Matrix{Any}(undef, 2, 1)
        tmp[1,1] = zero(I2)
        tmp[2,1] = adjoint(σ₋) * one(eltype(coefs))
        push!(data, SparseMPOTensor(tmp, T, pspace))
        positions = collect(cols[1]:row)
    else
        tmp = Matrix{Any}(undef, 1, 2)
        tmp[1, 1] = I2
        tmp[1, 2] = (-coefs[1]) * σ₊
        mpoj = SparseMPOTensor(tmp, T, pspace)
        data = [mpoj]
        positions = [cols[1]]

        for i in cols[1]+1:row-1
            tmp = Matrix{Any}(undef, 2, 2)
            tmp[1,1] = I2
            tmp[2,2] = m22JW
            tmp[2,1] = 0.
            pos = findfirst(x->x==i, cols)
            if isnothing(pos)
                tmp[1,2] = 0.
            else
                tmp[1,2] = (-coefs[pos]) * σ₊
            end
            push!(data, SparseMPOTensor(tmp, T, pspace))
            push!(positions, i)
        end
        tmp = Matrix{Any}(undef, 2, 2)
        tmp[1,1] = 0.
        tmp[2,2] = 0.
        tmp[1,2] = σ₊
        tmp[2,1] = adjoint(σ₋) * one(eltype(coefs))
        push!(data, SparseMPOTensor(tmp, T, pspace))
        push!(positions, row)
        for i in row+1 : cols[end]-1
            tmp = Matrix{Any}(undef, 2, 2)
            tmp[1,1] = I2
            tmp[2,2] = m22JW
            tmp[1,2] = 0.
            pos = findfirst(x->x==i, cols)
            if isnothing(pos)
                tmp[2,1] = 0.
            else
                tmp[2,1] = adjoint(σ₋) * coefs[pos]
            end
            push!(data, SparseMPOTensor(tmp, T, pspace))
            push!(positions, i)
        end
        tmp = Matrix{Any}(undef, 2, 1)
        tmp[1,1] = I2
        tmp[2,1] = adjoint(σ₋) * coefs[end]
        push!(data, SparseMPOTensor(tmp, T, pspace))
        push!(positions, cols[end])
    end

    mpo = MPO(MPOHamiltonian(data))
    return PartialMPO(mpo.data, positions)
end
