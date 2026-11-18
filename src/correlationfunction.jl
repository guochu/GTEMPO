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

# BCS bath：QuAPI 反常关联函数（cc/cn/nc/nn 四块有向指数系数；D7 仅 :z2 装配守卫在 bcs_partial_if）
# bcs_Δτ/Δt/Δm 自身返回 BCSCorrelationFunction（四味块 cc/cn/nc/nn 装配在 QuAPI 内完成）
correlationfunction(bath::AbstractBCSBath, lattice::ImagGrassmannLattice) =
    bcs_Δτ(bath.f; β=bath.β, μ=bath.μ, Δ=bath.Δ, N=lattice.Nτ, δτ=lattice.δτ)
correlationfunction(bath::AbstractBCSBath, lattice::RealGrassmannLattice) =
    bcs_Δt(bath.f; β=bath.β, μ=bath.μ, Δ=bath.Δ, N=lattice.Nt, t=lattice.t)
correlationfunction(bath::AbstractBCSBath, lattice::MixedGrassmannLattice) =
    bcs_Δm(bath.f; β=bath.β, μ=bath.μ, Δ=bath.Δ, Nτ=lattice.Nτ, t=lattice.t, Nt=lattice.Nt, δτ=lattice.δτ)

# ---------------------------------------------------------------------------
# coupling_norm：浴谱总强度（HybridizationKernel.coupling_norm 的单浴因子；
# 多带总量 = Σ_b coupling_norm(bath_b)。μ/β 不进谱权重——它们只进 Δ 装配的
# 核因子（spectrumshift(f, μ) 等），离散谱求和天然 μ 不变，bath 层直取谱）
# ---------------------------------------------------------------------------

"""
    coupling_norm(f) -> Real
    coupling_norm(bath::AbstractNormalBath) -> Real

浴谱总强度 α（= `hybriddynamics`/`sysdynamics_weyl2order` 的 `coupling_norm`
kwarg 在单带下的取值；多带总量 = `sum(coupling_norm, baths)`）：

- `DiracDelta`：α；
- `DiscreteSpectrum`：Σₙ αₙ（fs 为 δ-峰权重，vₙ = √αₙ = `spectrumcouplings`）；
- 连续 `AbstractBoundedFunction`：∫ f(ε)dε over [lb, ub]（quadgk 数值积分，
  容差同 IMB QuadGKTolerance）；
- δ×函数乘积（`DeltaMultF`/`DeltasMultF`，如热流的 ω 加权谱）：走各自
  quadgkwrapper 求和/取值规则（g ≡ 1）——同为 AbstractBoundedFunction，
  由 `quadgkwrapper(f)` 一并覆盖；
- `AbstractNormalBath`：直接取其谱——**与 μ、β 无关**；
- `AbstractBCSBath`：裸谱强度 ∫f dε（每带，`bath.spectrum` 直取——与 normal
  同款）。**Δ、μ、β 都不进**：Weyl 标量只数对角接触（āσaσ 槽）的宽度，
  Bogoliubov 正则性（u²+v²=1）保证反常通道（uv 块，ā↑ā↓/a↓a↑ 槽）对标量
  零贡献、总宽度 = 裸谱 normal 范数（推导与数值裁决见 `bcsweyl/README.md`：
  平凡系统接触缺口提取 cn*(Δ)/Σα ≡ 1（≤4e-4，Δ∈[0,2]）；generic H_S 下
  Z 收敛斜率 1.97–1.98）。一个 BCS 浴 = 一个带对（两带各 α）⟹ K 侧
  全带总量对该对按 2× 计入 `Σ_b coupling_norm(bath_b)`。

扩展点：自定义谱子类实现 `quadgkwrapper` 乘积路径或直接补 `coupling_norm`
方法。
"""
coupling_norm(f::DiracDelta) = f.α
coupling_norm(f::DiscreteSpectrum) = sum(spectrumvalues(f))
coupling_norm(f::AbstractBoundedFunction) = quadgkwrapper(f)
coupling_norm(bath::AbstractNormalBath) = coupling_norm(bath.spectrum)
coupling_norm(bath::AbstractBCSBath) = coupling_norm(bath.spectrum)
