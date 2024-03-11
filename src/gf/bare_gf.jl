###--------------imaginary time----------------
function Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...;
            band::Int=1, c1::Bool=false, c2::Bool=true, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = integrate(lattice, A, B..., alg=alg))
	pos1, pos2 = index(lattice, i, conj=c1, band=band), index(lattice, j, conj=c2, band=band)
	t = GTerm(pos1, pos2, coeff=1)
	A2 = _mult_A(t, A)
	return integrate(lattice, A2, B..., alg=alg)/Z
end
Gτ(lattice::ImagGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = Gτ(lattice, i, 1, A, B...; kwargs...)

function parallel_Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; 
            band::Int=1, c1::Bool=false, c2::Bool=true, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = parallel_integrate(lattice, A, B..., alg=alg) )
    pos1, pos2 = index(lattice, i, conj=c1, band=band), index(lattice, j, conj=c2, band=band)
    t = GTerm(pos1, pos2, coeff=1)
    A2 = _mult_A(t, A)
    return parallel_integrate(lattice, A2, B..., alg=alg)/Z
end
parallel_Gτ(lattice::ImagGrassmannLattice, i::Int, A::Vector{<:GrassmannMPS}, B::GrassmannMPS...; kwargs...) = parallel_Gτ(lattice, i, 1, A, B...; kwargs...)

function Gτ(lattice::ImagGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...;
            band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = integrate(lattice, A, B..., alg=alg))
    g = zeros(Float64, lattice.k)
    for i in 1:lattice.k-1
        g[i] = Gτ(lattice, i, A, B...; band=band, alg=alg, Z=Z)
    end
    g[end] = 1 - g[1]
    return g
end

function parallel_Gτ(lattice::ImagGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                     band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = parallel_integrate(lattice, A, B..., alg=alg))
    g = zeros(Float64, lattice.k)
    function _f(ist::Int, ifn::Int, r, lat, args...; kws...)
        for i in ist:ifn
            r[i] = Gτ(lat, i, args...; kws...)
        end
    end
    parallel_run(lattice.k-1, Threads.nthreads(), _f, g, lattice, A, B...; Z=Z, band=band, alg=alg)
    g[end] = 1 - g[1]
    return g
end

###--------------real time 1 order----------------
function Gt(lattice::RealGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
            b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Int=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg))
    pos1, pos2 = index(lattice, i, conj=c1, branch=b1, band=band), index(lattice, j, conj=c2, branch=b2, band=band)
    t = GTerm(pos1, pos2, coeff=1)
    A2 = _mult_A(t, A)
    g = integrate(lattice, A2, B..., alg=alg)/Z
    return g
end

###--------------mixed time 1 order----------------

"""
    Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...)

Mixed time Green's functions
"""
function Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
            b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Int=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg))
    pos1, pos2 = index(lattice, i, conj=c1, branch=b1, band=band), index(lattice, j, conj=c2, branch=b2, band=band)
    t = GTerm(pos1, pos2, coeff=1)
    A2 = _mult_A(t, A)
    g = integrate(lattice, A2, B..., alg=alg)/Z
    return g
end

function Gt(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; b1::Symbol, b2::Symbol, kwargs...)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    return Gm(lattice, i, j, A, B...; b1=b1, b2=b2, kwargs...)
end

function Gτ(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; c1::Bool=false, c2::Bool=true, kwargs...)
    return Gm(lattice, i, j, A, B...; b1=:τ, b2=:τ, c1=c1, c2=c2, kwargs...)
end
Gτ(lattice::MixedGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = Gτ(lattice, i, 1, A, B...; kwargs...)

function greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                 band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    @assert i >= j
    return Gt(lattice, i, j, A, B...; b1=:+, b2=:+, c1=false, c2=true, band=band, alg=alg, Z=Z)
end
greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = greater(lattice, i, 1, A, B...; kwargs...)
function lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                 band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    @assert i <= j
    return Gt(lattice, i, j, A, B...; b1=:-, b2=:+, c1=true, c2=false, band=band, alg=alg, Z=Z)
end
lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = lesser(lattice, 1, i, A, B...; kwargs...)


# other real-time observables
function occupation(lattice::RealGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) 
    return real(Gt(lattice, i, i, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end
occupation(lattice::RealGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; band::Int=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg)) = [occupation(lattice, i, A, B...; alg=alg, Z=Z, band=band) for i in 1:lattice.N]


# real-time first order
function electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                         alg::IntegrationAlgorithm=ExactIntegrate(),
                         Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1, max_range::Int=10000)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    for j in max(1, k-max_range):k-1
        curr += η⁺⁺[k, j] * Gt(lattice, k, j, A, B...; b1=:+, b2=:+, Z=Z, band=band, alg=alg)
        curr += η⁺⁻[k, j] * Gt(lattice, k, j, A, B...; b1=:+, b2=:-, Z=Z, band=band, alg=alg)
    end
    return 2 * curr / lattice.δt
end
electriccurrent(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; band::Int=1,
                alg::IntegrationAlgorithm=ExactIntegrate(), 
                Z::Number = integrate(lattice, A, B..., alg=alg)) = [electriccurrent(lattice, corr, k, A, B...; alg=alg, Z=Z, band=band) for k in 2:lattice.k]

"""
    electriccurrent_fast(A::GrassmannMPS, B::GrassmannMPS, lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int; kwargs...)

Calculate electric current by converting the current operator into an MPO
"""
function electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, k::Int, A::GrassmannMPS, B::GrassmannMPS...; band::Int=1, 
                            alg::IntegrationAlgorithm=ExactIntegrate(), 
                            Z::Number = integrate(lattice, A, B..., alg=alg))
    mpo = build_current_mpo(lattice, corr, k, band)
    A2 = mpo * A
    curr = integrate(lattice, A2, B..., alg=alg) / Z
    return 2 * curr / lattice.δt
end
electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr::RealCorrelationFunction, A::GrassmannMPS, B::GrassmannMPS...; band::Int=1,
                   alg::IntegrationAlgorithm=ExactIntegrate(),
                   Z::Real = integrate(lattice, A, B..., alg=alg)) = [electriccurrent_fast(lattice, corr, k, A, B...; alg=alg, Z=Z, band=band) for k in 2:lattice.k]


# real-time second order
function occupation(lattice::RealGrassmannLattice2Order, A::GrassmannMPS, B::GrassmannMPS...; kwargs...) 
    return real(Gt(lattice, lattice.k, lattice.k, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end

function electriccurrent(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::GrassmannMPS, B::GrassmannMPS...; 
                            alg::IntegrationAlgorithm=ExactIntegrate(), 
                            Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1, max_range::Int=10000)
    curr = complex(0.)
    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    @assert  lattice.N <= div(size(η⁺⁺,1)-1, 2) 
    k = lattice.k
    curr -= η⁺⁺[2*k-1, 1] * Gt(lattice, k, 1, A, B...; b1=:-, b2=:+, Z=Z, band=band, alg=alg)
    curr -= η⁺⁻[2*k-1, 1] * Gt(lattice, k, 1, A, B..., b1=:-, b2=:-, Z=Z, band=band, alg=alg)
    for j in 2:k-1
        curr -= (η⁺⁺[2*k-1, 2*j-2] + η⁺⁺[2*k-1, 2*j-1]) * Gt(lattice, k, j, A, B..., b1=:-, b2=:+, Z=Z, band=band, alg=alg)
        curr -= (η⁺⁻[2*k-1, 2*j-2] + η⁺⁻[2*k-1, 2*j-1])  * Gt(lattice, k, j, A, B..., b1=:-, b2=:-, Z=Z, band=band, alg=alg)
    end
    curr -= η⁺⁺[2*k-1, 2*k-2] * Gt(lattice, k, k, A, B..., b1=:-, b2=:+, Z=Z, band=band, alg=alg)
    curr -= η⁺⁻[2*k-1, 2*k-2] * Gt(lattice, k, k, A, B..., b1=:-, b2=:-, Z=Z, band=band, alg=alg)
    return 2 * curr / (0.5*lattice.δt)
end

function electriccurrent_fast(lattice::RealGrassmannLattice2Order, corr::RealCorrelationFunction, A::GrassmannMPS, B::GrassmannMPS...; 
                            alg::IntegrationAlgorithm=ExactIntegrate(), 
                            Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1)
    k = lattice.k
    mpo = build_current_mpo(lattice, corr, k, band)
    A2 = mpo * A
    curr = integrate(lattice, A2, B..., alg=alg) / Z

    η⁺⁺, η⁺⁻, η⁻⁺, η⁻⁻ = corr.G₊₊, corr.G₊₋, corr.G₋₊, corr.G₋₋
    curr -= η⁺⁺[2*k-1, 2*k-2] * Gt(lattice, k, k, A, B..., b1=:-, b2=:+, Z=Z, band=band, alg=alg)
    curr -= η⁺⁻[2*k-1, 2*k-2] * Gt(lattice, k, k, A, B..., b1=:-, b2=:-, Z=Z, band=band, alg=alg)
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

    virtual = isomorphism(Matrix{eltype(coefs)}, Rep[ℤ₂](1=>1), Rep[ℤ₂](1=>1))
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




