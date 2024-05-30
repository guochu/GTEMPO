Cτ(bath::AbstractFermionicBath; N::Int, δτ::Real=bath.β/N) = Cτ(bath.spectrum, β=bath.β, N=N, μ=bath.μ, δτ=δτ)
Cτ(f::SpectrumFunction; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) = Cτ(f, β, N, μ, δτ)


struct ImagCorrelationFunction{M<:AbstractMatrix{Float64}} <: AbstractCorrelationFunction
    data::M
end

Base.size(x::ImagCorrelationFunction, i::Int...) = size(x.data, i...)

function Base.show(io::IO, ::MIME"text/plain", A::ImagCorrelationFunction)
    print(io, "Imaginary Correlation Function [$(size(A.data, 1))]")
end

Base.:+(A::ImagCorrelationFunction, B::ImagCorrelationFunction) = ImagCorrelationFunction(A.data + B.data)
index(x::ImagCorrelationFunction, i::Int, j::Int) = x.data[i, j]


"""
    Cτ(f, β::Real, N::Int)

f is the spectrum function
"""
function Cτ(f0::SpectrumFunction, β::Real, N::Int, μ::Real, δτ::Real=β / N)
    f′, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    f(ϵ) = f′(ϵ + μ)
    lb -= μ
    ub -= μ
    # δτ = β / N
    
    g₁(ϵ) = _g₁(β, 0., ϵ)
    g₂(ϵ) = _g₂(β, 0., ϵ)
    fⱼₖ(Δk::Int, ε::Float64) = _fⱼₖ_i(f, Δk, ε, δτ)
    fₖⱼ(Δk::Int, ε::Float64) = _fₖⱼ_i(f, Δk, ε, δτ)
    fⱼⱼ(ε::Float64) = _fⱼⱼ_i(f, ε, δτ)
    fₖₖ(ε::Float64) = _fₖₖ_i(f, ε, δτ)

    # j >= k
    L = N
    ηⱼₖ = zeros(Float64, L)
    ηⱼₖ[1] = quadgk(ε -> g₁(ε)*fⱼⱼ(ε), lb, ub)[1]
    for k = 1:L-1
        ηⱼₖ[k+1] = quadgk(ε -> g₁(ε)*fⱼₖ(k,ε), lb, ub)[1]
    end

    ηₖⱼ = zeros(Float64, L)
    ηₖⱼ[1] = quadgk(ε -> -g₂(ε)*fₖₖ(ε), lb, ub)[1]
    for k = 1:L-1
        ηₖⱼ[k+1] = quadgk(ε -> -g₂(ε)*fₖⱼ(k,ε), lb, ub)[1]
    end
    ImagCorrelationFunction(CorrelationMatrix{Float64}(ηⱼₖ, ηₖⱼ))
end

# from now on j > k
function _fⱼₖ_i(f, Δk::Int, ε::Float64, δτ)
    if (abs(ε) > tol)
        -2f(ε)/ε^2*exp(-Δk*δτ*ε)*(1-cosh(δτ*ε))
    else
        f(ε)*exp(-Δk*δτ*ε)*δτ^2
    end
end

function _fₖⱼ_i(f, Δk::Int, ε::Float64, δτ)
    if (abs(ε) > tol)
        -2f(ε)/ε^2*exp(Δk*δτ*ε)*(1-cosh(δτ*ε))
    else
        f(ε)*exp(Δk*δτ*ε)*δτ^2
    end
end

function _fⱼⱼ_i(f, ε::Float64, δτ)
    if (abs(ε) > tol)
        f(ε)/ε^2*(exp(-δτ*ε)-(1-δτ*ε))
    else
        0.5f(ε)*δτ^2
    end
end

function _fₖₖ_i(f, ε::Float64, δτ)
    if (abs(ε) > tol)
        f(ε)/ε^2*(exp(δτ*ε)-(1+δτ*ε))
    else
        0.5*f(ε)*δτ^2/2
    end
end

function Δiw_to_Cτ(Δiw::AbstractVector{<:Number}; β::Real, N::Int)
    iseven(length(Δiw)) || throw(ArgumentError("even number of frequencies expected"))
    δτ = β/N
    nmax = div(length(Δiw), 2) - 1

    ηⱼₖ = zeros(ComplexF64, N)
    ηₖⱼ = zeros(ComplexF64, N)

    for (i, n) in enumerate(-nmax:nmax+1)
        ωₙ = (2n-1)*π/β
        a = -(2/(β*ωₙ^2)) * Δiw[i] * (1 - cos(ωₙ * δτ))
        ηⱼₖ[1] += a / 2
        ηₖⱼ[1] += a / 2
        for k in 1:N-1
            ηⱼₖ[k+1] += a * exp(-im*ωₙ*k*δτ)
            ηₖⱼ[k+1] += a * exp(im*ωₙ*k*δτ)
        end
    end
    for i in 1:N
        (abs(imag(ηⱼₖ[i])) < 1.0e-8) || error("imag part of ηⱼₖ[$i] is nonzero")
        (abs(imag(ηₖⱼ[i])) < 1.0e-8) || error("imag part of ηₖⱼ[$i] is nonzero")
    end

    ImagCorrelationFunction(CorrelationMatrix{Float64}(real(ηⱼₖ), real(ηₖⱼ)))
end