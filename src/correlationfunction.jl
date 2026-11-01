# wrapper for QuAPI package

"""
    correlationfunction(bath::AbstractBath, lattice::AbstractGrassmannLattice)

Compute the discrete hybridization functions using QuAPI
"""
function correlationfunction(bath::AbstractBath, lattice::ImagGrassmannLattice)
    # @assert lattice.β == bath.β
    isapprox(lattice.β, bath.β; rtol=1.0e-12, atol=1.0e-14) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δτ(bath, N=lattice.N, δτ=lattice.δτ)
end 
correlationfunction(bath::AbstractBath, lattice::RealGrassmannLattice) = Δt(bath, N=lattice.N, t=lattice.t) 
function correlationfunction(bath::AbstractBath, lattice::MixedGrassmannLattice)
    isapprox(lattice.β, bath.β; rtol=1.0e-12, atol=1.0e-14) || @warn "lattice.β=$(lattice.β), but bath.β=$(bath.β)"
    Δm(bath, Nτ=lattice.Nτ, t=lattice.t, Nt=lattice.Nt)
end  

# ImagCorrelationFunction{T} 的参数 T 是 CorrelationMatrix 的元素类型（QuAPI 2026-10-30 重构后）
TK.scalartype(::Type{<:ImagCorrelationFunction{T}}) where {T} = T
TK.scalartype(::Type{<:RealCorrelationFunction}) = ComplexF64
TK.scalartype(::Type{<:MixedCorrelationFunction}) = ComplexF64
TK.scalartype(::Type{<:BCSCorrelationFunction{M}}) where M = scalartype(M)
# ---------------------------------------------------------------------------
# QuAPI 统一实时约定注记（2026-10-30）：两支 Keldysh 关联函数
#（RealCorrelationFunction，含 BCS 包装在实时轮廓上的 uu/ud/du）把跨支块的
# Keldysh 相对负号嵌在块值里（(+,−) = −(c, b)、(−,+) = −(a, d)；同支块不变）。
# GTEMPO 的 IF 指数与可观测量公式消费的是无符号的物理块，因此**跨支读数处
# 要补一个负号**——这个约定由各调用点就地处理（partialif 的系数读取、
# ttiif/realtime 的分支装配、currents 的 η⁺⁻），判定条件统一用
# `branches(corr) == (:+, :-) && b1 !== b2`（branches 对 BCS 委托到 uu，
# 一个判断同时覆盖普通与 BCS）。虚时 τ 环与混合轮廓的块值即物理值，不动；
# 玻色（Fock 路径）的历史存储本就带符号，同样不动（见 fockcontourindex.jl）。
# ---------------------------------------------------------------------------
