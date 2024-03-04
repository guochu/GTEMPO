abstract type AbstractCorrelationFunction end


struct CorrelationMatrix{T <: Number} <: AbstractMatrix{T}
    ηⱼₖ::Vector{T}  # j >= k
    ηₖⱼ::Vector{T}    
end

Base.size(x::CorrelationMatrix, i::Int) = ifelse(i <= 2, length(x.ηⱼₖ), 1)
Base.size(x::CorrelationMatrix) = (length(x.ηⱼₖ), length(x.ηⱼₖ))

Base.similar(a::CorrelationMatrix, dims::Tuple{Int, Int}) = Matrix{eltype(a)}(undef, dims)
Base.similar(a::CorrelationMatrix, T::Type, dims::Tuple{Int, Int}) = Matrix{T}(undef, dims)

function Base.getindex(A::CorrelationMatrix, i::Int, j::Int)
    L = size(A, 1)
    ((1 <= i <= L) && (1 <= j <= L)) || throw(BoundsError())
    if (i > j)
        A.ηⱼₖ[i-j+1]
    elseif (i < j)
        A.ηₖⱼ[j-i+1]
    else
        A.ηⱼₖ[1] + A.ηₖⱼ[1]
    end
end

function Base.:+(x::CorrelationMatrix, y::CorrelationMatrix)
    (size(x) == size(y)) || throw(DimensionMismatch("correlation matrix size mismatch"))
    return CorrelationMatrix(x.ηⱼₖ + y.ηⱼₖ, x.ηₖⱼ + y.ηₖⱼ)
end
Base.:(-)(x::CorrelationMatrix) = CorrelationMatrix(-x.ηⱼₖ, -x.ηₖⱼ)
Base.transpose(x::CorrelationMatrix) = CorrelationMatrix(x.ηₖⱼ, x.ηⱼₖ)

include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")
include("mixedtime.jl")


"""
    correlationfunction(bath::AbstractFermionicBath, lattice::AbstractGrassmannLattice)

Compute the discrete correlation functions (QUAPI)
"""
function correlationfunction(bath::AbstractFermionicBath, lattice::ImagGrassmannLattice1Order)
    # @assert lattice.β == bath.β
    (lattice.β == bath.β) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δτ(bath, N=lattice.N, δτ=lattice.δτ)
end 
correlationfunction(bath::AbstractFermionicBath, lattice::RealGrassmannLattice1Order) = Δt(bath, N=lattice.N, t=lattice.t) 
correlationfunction(bath::AbstractFermionicBath, lattice::RealGrassmannLattice2Order) = Δt(bath, N=2*lattice.N, t=lattice.t)