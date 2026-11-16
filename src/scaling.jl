# scaling.jl：缩放标量（后端直用）与有界缓存。
#
# 缩放标量 = **GrassmannMPS.ScaledScalar 直接使用**（2026-11-16：旧 ScaledComplex
# 包装退役——GM.ScaledScalar 为其严格超集：参数化实/复相位（实链结果保持实型）、
# 同构造归一化、value ≡ 旧 unscale（同 OverflowError 守卫）、完整 scaled 算术
# （×/÷/+/−/conj/abs，log 域精确）。公共 API：ScaledScalar/value 由主文件
# using+export 透传；unscale = value 的历史名别名；scaled_ratio = 精确比值物化。
const unscale = GM.value
scaled_ratio(a::GM.ScaledScalar, b::GM.ScaledScalar) = GM.value(a / b)

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
