# scaling.jl：缩放标量与有界缓存（U1GTEMPO2 scaling.jl 原样移植）

"A complex scalar represented by log(abs(value)) and a unit complex phase."
struct ScaledComplex
    logabs::Float64
    phase::ComplexF64
    function ScaledComplex(logabs::Real, phase::Number)
        (isfinite(logabs) || logabs == -Inf) && isfinite(phase) || throw(ArgumentError("nonfinite scaled scalar"))
        if iszero(phase) || logabs == -Inf
            new(-Inf, 0im)
        else
            new(Float64(logabs) + log(abs(phase)), ComplexF64(phase / abs(phase)))
        end
    end
end
ScaledComplex(z::Number) = ScaledComplex(0.0, z)
# 公共 API 返回 ScaledComplex；内部与 GrassmannMPS ScaledScalar 双向互转
ScaledComplex(z::GM.ScaledScalar) = ScaledComplex(z.logabs, ComplexF64(z.phase))
Base.iszero(z::ScaledComplex) = z.logabs == -Inf
function unscale(z::ScaledComplex)
    iszero(z) && return 0.0 + 0im
    log(nextfloat(0.0)) <= z.logabs <= log(floatmax(Float64)) ||
        throw(OverflowError("complex scalar outside Float64 range; use the scaled_* channels"))
    z.phase * exp(z.logabs)
end
function scaled_ratio(a::ScaledComplex, b::ScaledComplex)
    iszero(b) && throw(ArgumentError("zero denominator"))
    iszero(a) && return 0.0 + 0im
    unscale(ScaledComplex(a.logabs - b.logabs, a.phase / b.phase))
end

"Bounded FIFO cache owned by one preparation or solution, never process-global."
mutable struct BoundedCache
    values::Dict{Any,Any}
    order::Vector{Any}
    bytes::Int
    maxentries::Int
    maxbytes::Int
    hits::Int
    misses::Int
end
function BoundedCache(; maxentries=256, maxbytes=64 * 1024^2)
    maxentries >= 0 && maxbytes >= 0 || throw(ArgumentError("negative cache bound"))
    BoundedCache(Dict{Any,Any}(), Any[], 0, maxentries, maxbytes, 0, 0)
end
Base.@noinline function cached!(f, cache::BoundedCache, key)
    Base.@nospecialize key
    if haskey(cache.values, key)
        cache.hits += 1; return cache.values[key][1]
    end
    cache.misses += 1; value = f(); bytes = Base.summarysize((key, value))
    if cache.maxentries > 0 && bytes <= cache.maxbytes
        while !isempty(cache.order) && (length(cache.order) >= cache.maxentries || cache.bytes + bytes > cache.maxbytes)
            old = popfirst!(cache.order); cache.bytes -= cache.values[old][2]; delete!(cache.values, old)
        end
        cache.values[key] = (value, bytes); push!(cache.order, key); cache.bytes += bytes
    end
    value
end
function Base.empty!(c::BoundedCache)
    empty!(c.values); empty!(c.order); c.bytes = 0; c.hits = 0; c.misses = 0; c
end
