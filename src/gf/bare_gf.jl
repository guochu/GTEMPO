###--------------imaginary time----------------

function gf(lattice::AbstractGrassmannLattice, a::NTuple{N, ContourIndex}, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg)) where {N}
    pos = map(x->lattice[x], a)
    t = GTerm(pos, coeff=1)
    A2 = _mult_A(t, A)
    return integrate(lattice, A2, B..., alg=alg)/Z    
end

"""
    gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Calculate the Green's function ⟨x y⟩, where x, y are GVs 
specified by contour indices a and b

If A is type Vector{GrassmannMPS}, it means A is a sum of GMPSs
The multiplication of A and B... is assumed
alg: the algorithm to perform the integration
Z: the partition function values
"""
function gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    pos1, pos2 = lattice[a], lattice[b]
    t = GTerm(pos1, pos2, coeff=1)
    A2 = _mult_A(t, A)
    return integrate(lattice, A2, B..., alg=alg)/Z    
end

"""
    contour_ordered_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
                alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

similar to gf, but return the contour ordered Green's function
"""
function contour_ordered_gf(lattice::AbstractGrassmannLattice, a::ContourIndex, b::ContourIndex, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
                            alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = integrate(lattice, A, B..., alg=alg)) 
    ((!a.conj) && (b.conj)) || throw(ArgumentError("conj(a)=false and conj(b)=true should be satisfied"))
    return (a < b) ? -gf(lattice, b, a, A, B...; alg=alg, Z=Z) : gf(lattice, a, b, A, B...; alg=alg, Z=Z) 
end

"""
    Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
        band::Int=1, c1::Bool=false, c2::Bool=true, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the Matsubara Green's function ⟨aᵢ bⱼ⟩ on the imaginary-time axis
band: the band to calculate the Green's function
c1: whether to take the conjugate of the first GV a
c2: whether to take the conjugate of the second GV b
The other keywords are the same as gf
"""
function Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Int=1, c1::Bool=false, c2::Bool=true, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Real = integrate(lattice, A, B..., alg=alg))
    a, b = ContourIndex(i, conj=c1, branch=:τ, band=band), ContourIndex(j, conj=c2, branch=:τ, band=band)
    return gf(lattice, a, b, A, B...; alg=alg, Z=Z)
end
Gτ(lattice::ImagGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = Gτ(lattice, i, 1, A, B...; kwargs...)

###--------------real time 1 order----------------

"""
    Gt(lattice::RealGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
        b1, b2, band::Int=1, c1::Bool=true, c2::Bool=false, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the real-time Green's function ⟨aᵢ bⱼ⟩ on the Keldysh contour
band: the band to calculate the Green's function
c1: whether to take the conjugate of the first GV a
c2: whether to take the conjugate of the second GV b
b1: the branch of the first GV (can be :+ or :-)
b2: the branch of the second GV (can be :+ or :-)
The other keywords are the same as gf
"""
function Gt(lattice::RealGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS}; 
            b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Int=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg))
    a, b = ContourIndex(i, conj=c1, branch=b1, band=band), ContourIndex(j, conj=c2, branch=b2, band=band)
    return gf(lattice, a, b, A, B...; alg=alg, Z=Z)
end

###--------------mixed time 1 order----------------

"""
    Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
        b1, b2, band::Int=1, c1::Bool=true, c2::Bool=false, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the Green's function ⟨aᵢ bⱼ⟩ on the L-shaped Kadanoff-Baym contour
band: the band to calculate the Green's function
c1: whether to take the conjugate of the first GV a
c2: whether to take the conjugate of the second GV b
b1: the branch of the first GV (can be :+, :- or :τ)
b2: the branch of the second GV (can be :+, :- or :τ)
The other keywords are the same as gf
"""
function Gm(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
            b1::Symbol, b2::Symbol, c1::Bool=true, c2::Bool=false, band::Int=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg))
    
    a, b = ContourIndex(i, conj=c1, branch=b1, band=band), ContourIndex(j, conj=c2, branch=b2, band=band)
    return gf(lattice, a, b, A, B...; alg=alg, Z=Z)
end

########************************#######

"""
    parallel_Gτ(lattice::ImagGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
        band::Int=1, c1::Bool=false, c2::Bool=true, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

The same as Gτ, but use multi-threading
"""
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

function Gt(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; b1::Symbol, b2::Symbol, kwargs...)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    return Gm(lattice, i, j, A, B...; b1=b1, b2=b2, kwargs...)
end

function Gτ(lattice::MixedGrassmannLattice, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; c1::Bool=false, c2::Bool=true, kwargs...)
    return Gm(lattice, i, j, A, B...; b1=:τ, b2=:τ, c1=c1, c2=c2, kwargs...)
end
Gτ(lattice::MixedGrassmannLattice, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = Gτ(lattice, i, 1, A, B...; kwargs...)

"""
    greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Int=1, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the greater Green's function ⟨aᵢ bⱼ⟩ on the Keldysh or the Kadanoff-Baym contour
"""
function greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                 band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    @assert i >= j
    return Gt(lattice, i, j, A, B...; b1=:+, b2=:+, c1=false, c2=true, band=band, alg=alg, Z=Z)
end
greater(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = greater(lattice, i, 1, A, B...; kwargs...)

"""
    lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Int=1, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the lesser Green's function ⟨aᵢ bⱼ⟩ on the Keldysh or the Kadanoff-Baym contour
"""
function lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, j::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                 band::Int=1, alg::IntegrationAlgorithm=ExactIntegrate(), Z::Number = integrate(lattice, A, B..., alg=alg))
    @assert i <= j
    return Gt(lattice, i, j, A, B...; b1=:-, b2=:+, c1=true, c2=false, band=band, alg=alg, Z=Z)
end
lesser(lattice::Union{RealGrassmannLattice, MixedGrassmannLattice}, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) = lesser(lattice, 1, i, A, B...; kwargs...)


# other real-time observables

"""
    occupation(lattice::RealGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::Vararg{GrassmannMPS};
            band::Int=1, alg=ExactIntegrate(), Z = integrate(lattice, A, B..., alg=alg))

Return the occupation at time step i on the Keldysh contour
"""
function occupation(lattice::RealGrassmannLattice1Order, i::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; kwargs...) 
    return real(Gt(lattice, i, i, A, B...; c1=false, c2=true, b1=:+, b2=:-, kwargs...))
end
occupation(lattice::RealGrassmannLattice1Order, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; band::Int=1, 
            alg::IntegrationAlgorithm=ExactIntegrate(), 
            Z::Number = integrate(lattice, A, B..., alg=alg)) = [occupation(lattice, i, A, B...; alg=alg, Z=Z, band=band) for i in 1:lattice.N]


# real-time first order

"""
    electriccurrent(lattice::RealGrassmannLattice1Order, corr, k::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                         alg::IntegrationAlgorithm=ExactIntegrate(),
                         Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1, max_range::Int=10000)

Return the electric current at time step k, which is calculated as the sum
of a series of Green's functions
"""
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
    electriccurrent_fast(lattice::RealGrassmannLattice1Order, corr, k::Int, A::Union{GrassmannMPS, Vector}, B::GrassmannMPS...; 
                         alg::IntegrationAlgorithm=ExactIntegrate(),
                         Z::Number = integrate(lattice, A, B..., alg=alg), band::Int=1, max_range::Int=10000)

The same as electriccurrent, but use a very efficient algorithm which directly 
builds the current operator as an MPO of bond dimension 2
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




