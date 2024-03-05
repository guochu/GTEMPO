"""
	struct MixedCorrelationFunction

Correlations on the L-shaped contour,
there are 9 of them due to the double integral
"""
struct MixedCorrelationFunction <: AbstractCorrelationFunction
    # real time
    ηⱼₖ::Vector{ComplexF64}
    ηₖⱼ::Vector{ComplexF64}
    # imag time
    ξⱼₖ::Vector{ComplexF64}
    ξₖⱼ::Vector{ComplexF64}
    # mix time
    ζⱼₖ::Matrix{ComplexF64}
    ζₖⱼ::Matrix{ComplexF64}
end

function Base.show(io::IO, ::MIME"text/plain", A::MixedCorrelationFunction)
    print(io, "Mixed Correlation Function [$(isize(A))+$(rsize(A))]")
end

Base.:+(A::MixedCorrelationFunction, B::MixedCorrelationFunction) = MixedCorrelationFunction(A.ηⱼₖ + B.ηⱼₖ, A.ηₖⱼ + B.ηₖⱼ, A.ξⱼₖ + B.ξⱼₖ, A.ξₖⱼ + B.ξₖⱼ, A.ζⱼₖ + B.ζⱼₖ, A.ζₖⱼ + B.ζₖⱼ)
# branch(x::MixedCorrelationFunction, f1::Symbol, f2::Symbol) = ifelse(f1, ifelse(f2, x.G₊₊, x.G₊₋), ifelse(f2, x.G₋₊, x.G₋₋))

isize(x::MixedCorrelationFunction) = length(x.ξⱼₖ)
rsize(x::MixedCorrelationFunction) = length(x.ηⱼₖ)

Δm(bath::AbstractFermionicBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = Δm(bath.spectrum, β=bath.β, μ=bath.μ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)
function Δm(f0::SpectrumFunction; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
	f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    δt = t / Nt
    g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
    # real time
    fⱼₖ(Δk, ε) = _fⱼₖ_r(f, Δk, ε, δt)
    fⱼⱼ(ε) = _fⱼⱼ_r(f, ε, δt)
    fₖₖ(ε) = _fₖₖ_r(f, ε, δt)
    # imaginary time
    hⱼₖ(Δk::Int, ε::Float64) = _fⱼₖ_i(f, Δk, ε, δτ)
    hⱼⱼ(ε::Float64) = _fⱼⱼ_i(f, ε, δτ)
    hₖₖ(ε::Float64) = _fₖₖ_i(f, ε, δτ)
    # mixed time
    lⱼₖ(j::Int, k::Int, ε::Float64) = _lⱼₖ(f, j, k, ε, δt, δτ)
    lₖⱼ(k::Int, j::Int, ε::Float64) = _lₖⱼ(f, k, j, ε, δt, δτ)
    N, M = Nt, Nτ

    # real time
    ηⱼₖ = zeros(ComplexF64, Nt+1)
    ηₖⱼ = zeros(ComplexF64, Nt+1)

    ηⱼₖ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε), lb, ub)[1]
    for i = 1:Nt
        ηⱼₖ[i+1] = quadgk(ε -> -g₁(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> g₂(ε)*fₖₖ(ε), lb, ub)[1]
    for i = 1:Nt
        ηₖⱼ[i+1] = quadgk(ε -> g₂(ε)*fⱼₖ(-i,ε), lb, ub)[1]
    end

    # imag time
    ξⱼₖ = zeros(ComplexF64, M) # j >= k
    ξₖⱼ = zeros(Float64, M)

    ξⱼₖ[1] = quadgk(ε -> -g₁(ε)*hⱼⱼ(ε), lb, ub)[1]
    for k = 2:M
        ξⱼₖ[k] = quadgk(ε -> -g₁(ε)*hⱼₖ(k-1,ε), lb, ub)[1]
    end
    
    ξₖⱼ[1] = quadgk(ε -> g₂(ε)*hₖₖ(ε), lb, ub)[1]
    for k = 2:M
        ξₖⱼ[k] = quadgk(ε -> g₂(ε)*hⱼₖ(1-k,ε), lb, ub)[1]
    end

    # mix time
    ζⱼₖ = zeros(ComplexF64, M, N+1)
    ζₖⱼ = zeros(ComplexF64, N+1, M)

    for j = 1:M, k = 1:N+1
        ζⱼₖ[j,k] = quadgk(ε -> im*g₁(ε)*lⱼₖ(j-1,k-1,ε), lb, ub)[1]
        ζₖⱼ[k,j] = quadgk(ε -> -im*g₂(ε)*lₖⱼ(k-1,j-1,ε), lb, ub)[1]
    end
    MixedCorrelationFunction(ηⱼₖ,ηₖⱼ,ξⱼₖ,ξₖⱼ,ζⱼₖ,ζₖⱼ)
end

function index(A::MixedCorrelationFunction, i::Int, j::Int; b1::Symbol, b2::Symbol)
	b₁, b₂ = b1, b2
	@boundscheck begin
        (b₁ in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
        (b₂ in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))

		if b₁ == :τ
			(1 <= i <= isize(A)) || throw(BoundsError(1:isize(A), i))
		else
			(1 <= i <= rsize(A)) || throw(BoundsError(1:rsize(A), i))
		end
		if b₂ == :τ
			(1 <= j <= isize(A)) || throw(BoundsError(1:isize(A), j))
		else
			(1 <= j <= rsize(A)) || throw(BoundsError(1:rsize(A), j))
		end
	end
    # here b₁ and b₂ are branches
    if (b₁ == :+ && b₂ == :+) # G¹¹=G⁺⁺
        if (i == j)
        	A.ηⱼₖ[1] + A.ηₖⱼ[1]
        elseif (i > j)
            A.ηⱼₖ[i-j+1]
        else
            A.ηₖⱼ[j-i+1]
        end
    elseif (b₁ == :+ && b₂ == :-) # G¹²=G⁺⁻
        if (i == j)
            -2*A.ηₖⱼ[1]
        elseif (i > j)
            -A.ηₖⱼ[i-j+1]'
        else
            -A.ηₖⱼ[j-i+1]
        end
    elseif (b₁ == :+ && b₂ == :τ) # G¹³
        -im*A.ζₖⱼ[i,j]
    elseif (b₁ == :- && b₂ == :+) # G²¹=G⁻⁺
        if (i == j)
            -2*A.ηⱼₖ[1]
        elseif (i > j)
            -A.ηⱼₖ[i-j+1]
        else
            -A.ηⱼₖ[j-i+1]'
        end
    elseif (b₁ == :- && b₂ == :-) # G²²=G⁻⁻
        if (i == j)
            A.ηⱼₖ[1]+A.ηₖⱼ[1]
        elseif (i > j)
            A.ηₖⱼ[i-j+1]
        else
            A.ηⱼₖ[j-i+1]
        end
    elseif (b₁ == :- && b₂ == :τ) # G²³
        im*A.ζₖⱼ[i,j]
    elseif (b₁ == :τ && b₂ == :+) # G³¹
        -im*A.ζⱼₖ[i,j]
    elseif (b₁ == :τ && b₂ == :-) # G³²
        im*A.ζⱼₖ[i,j]
    elseif (b₁ == :τ && b₂ == :τ) # G³³=G(τ)
        if (i == j)
            -A.ξⱼₖ[1]+A.ξₖⱼ[1]
        elseif (i > j)
            -A.ξⱼₖ[i-j+1]
        else
            -A.ξₖⱼ[j-i+1]
        end
    end
end


# mixed part lⱼₖ for G31 and lₖⱼ for G13
function _lⱼₖ(f, j::Int, k::Int, ε::Float64, δt, δτ)
    if (abs(ε) > tol)
        im*f(ε)/ε^2*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)
    else
        f(ε)*exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt
    end
end

function _lₖⱼ(f, k::Int, j::Int, ε::Float64, δt, δτ)
    if (abs(ε) > tol)
        -im*f(ε)/ε^2*exp(ε*j*δτ)*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)
    else
        -f(ε)*exp(ε*j*δτ)*exp(-im*ε*k*δt)*δτ*δt
    end
end

