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

include("util.jl")
include("imaginarytime.jl")
include("realtime.jl")