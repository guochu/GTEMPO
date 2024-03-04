Gt(bath::AbstractFermionicBath; N::Int, t::Real) = Gt(bath.spectrum, β=bath.β, N=N, t=t, μ=bath.μ)
Gt(f::SpectrumFunction; β::Real, N::Int, t::Real, μ::Real=0) = Gt(f, β, N, t/N, μ)

struct RealCorrelationFunction{A<:AbstractMatrix{ComplexF64}, B<:AbstractMatrix{ComplexF64}, C<:AbstractMatrix{ComplexF64}, D<:AbstractMatrix{ComplexF64}} <: AbstractCorrelationFunction
    G₊₊::A
    G₊₋::B
    G₋₊::C
    G₋₋::D
end

function Base.show(io::IO, ::MIME"text/plain", A::RealCorrelationFunction)
    print(io, "Real Correlation Function [$(size(A.G₊₊, 1))]")
end
index(x::RealCorrelationFunction, i::Int, j::Int; f1::Bool, f2::Bool) = branch(x, f1, f2)[i, j]

Base.:+(A::RealCorrelationFunction, B::RealCorrelationFunction) = RealCorrelationFunction(A.G₊₊ + B.G₊₊, A.G₊₋ + B.G₊₋, A.G₋₊ + B.G₋₊, A.G₋₋ + B.G₋₋)
branch(x::RealCorrelationFunction, f1::Bool, f2::Bool) = ifelse(f1, ifelse(f2, x.G₊₊, x.G₊₋), ifelse(f2, x.G₋₊, x.G₋₋))
function branch(x::RealCorrelationFunction, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    if b1 == :+
        if b2 == :+
            return x.G₊₊
        else
            return x.G₊₋
        end
    else
        if b2 == :+
            return x.G₋₊
        else
            return x.G₋₋
        end
    end
end
branch(x::RealCorrelationFunction; b1::Symbol, b2::Symbol) = branch(x, b1, b2)

function Gt(f0::SpectrumFunction, β::Real, N::Int, δt::Real, μ::Real)
    f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
    fⱼₖ(Δk, ε) = _fⱼₖ_r(f, Δk, ε, δt)
    fⱼⱼ(ε) = _fⱼⱼ_r(f, ε, δt)

    ### G₊₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε), lb, ub)[1] 
    for i = 1:N
        ηⱼₖ[i+1] = quadgk(ε -> -g₁(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> g₂(ε)*fⱼⱼ(ε)', lb, ub)[1]
    for i = 1:N
        ηₖⱼ[i+1] = quadgk(ε -> g₂(ε)*fⱼₖ(i,ε)', lb, ub)[1]
    end
    G₊₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)   

    ### G₊₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgk(ε -> g₂(ε)*fⱼⱼ(ε), lb, ub)[1] 
    for i = 1:N
        ηⱼₖ[i+1] = quadgk(ε -> g₂(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> g₂(ε)*fⱼⱼ(ε)', lb, ub)[1]
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₊₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε), lb, ub)[1] 
    for i = 1:N
        ηⱼₖ[i+1] = quadgk(ε -> -g₁(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε)', lb, ub)[1]
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₋₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgk(ε -> g₂(ε)*fⱼⱼ(ε), lb, ub)[1] 
    for i = 1:N
        ηⱼₖ[i+1] = quadgk(ε -> g₂(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε)', lb, ub)[1]
    for i = 1:N
        ηₖⱼ[i+1] = quadgk(ε -> -g₁(ε)*fⱼₖ(i,ε)', lb, ub)[1]
    end
    G₋₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)
    return RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
end


# assume j > k
function _fⱼₖ_r(f, Δk::Int, ε::Float64, δt)
    if (abs(ε) > tol)
        2*f(ε)/ε^2*exp(-im*ε*Δk*δt)*(1-cos(ε*δt))
    else
        f(ε)*exp(-im*ε*Δk*δt)*δt^2
    end
end

function _fⱼⱼ_r(f, ε::Float64, δt)
    if (abs(ε) > tol)
        f(ε)/ε^2*((1-im*ε*δt)-exp(-im*ε*δt))
    else
        0.5*f(ε)*δt^2
    end
end

