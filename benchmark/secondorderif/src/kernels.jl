# ============================================================================
# kernels.jl — Weyl N×N 浴影响核（单元级）与耦合范数
#
# IF_raw = Z_bath · exp(Φ), Φ = Σ_{iα,jγ} ā_{iα} Δ_{iα,jγ} a_{jγ}，
# Δ = 原浴二点核在单元上的精确双重积分：跨单元为矩形区域，同单元保留两个
# 有序三角区域（推导：benchmark/weyl-secondorder.html）。
#
# 实时核（关键）：QuAPI `Δt` 的原始值是**顶点约定**（_vs_rt_acc 映射下的
# 节点读数），不是 Weyl 全宽单元的积分值。这里对单模 Dirac-δ 浴直接从
# **精确浴迹多项式**提取单元级核：单元演化算符
#     E_j^± = e^{∓iδt·H_b}·[I + G ± G²/2]，G = ∓iδt(ā F' + F'† a)
# （F' 为海森堡共轭后的浴湮灭矩阵，G³ = 0 使级数精确闭合），按轮廓时序
# (+ 支 1..N，再 − 支 N..1) 乘积并取浴迹，其线性单项式系数恰为
# exp(Φ) 的 Φ 核。交叉验证：τ 核直接用 QuAPI `Δτ`（本身按单元给出，
# 已与 ED 验证一致）。
#
# 块值沿用 QuAPI 的存储约定（跨支块内嵌 Keldysh 相对负号），因此可以逐块
# 交给 GTEMPO 的 partialif_hybrid/hybriddynamics 消费（它在跨支读数处补
# 负号）；同时除以浴迹常数（QuAPI 核的归一化约定，归一化观测量下抵消）。
# ============================================================================
using GTEMPO: Δt, Δτ, Δm, fermionicbath
using ImpurityModelBase: DiscreteSpectrum, spectrumvalues

"""
    HybridizationKernel

单元级杂化核：`blocks[(b1,b2)]` 的 (i,j) 元 = ā_i^{b1} 与 a_j^{b2} 之间的
IF 指数系数（含 QuAPI 存储约定的跨支 Keldysh 负号），单元下标 1..N。
`coupling_norm = Σ_k |V_k|²` 是 Weyl 标量补偿所需的单带总耦合权重
（U1GTEMPO2：order=2 必须提供 coupling_norm；有限采样核无法可靠推断）。
"""
struct HybridizationKernel
    contour::AbstractWeylContour
    blocks::Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}
    coupling_norm::Float64
end

coupling_norm(k::HybridizationKernel) = k.coupling_norm

# DiracDelta 谱的积分权重 α（精确）；一般有界谱数值积分。单模基准只需前者。
α_weight(spec::DiracDelta) = spec.α
# 离散（多 δ）谱：总耦合 = Σ αₙ（quadgk 对 ±Inf 上的离散谱无意义）
α_weight(spec::DiscreteSpectrum) = sum(spectrumvalues(spec))
α_weight(spec) = quadgk(spec, lowerbound(spec), upperbound(spec))[1]
α_weight(bath::AbstractBath) = α_weight(bath.spectrum)

"""
    hybridization(bath, c::AbstractWeylContour; band=1) -> HybridizationKernel

单元级浴核（三轮廓统一 QuAPI 谱路径）：
- 虚时：QuAPI `Δτ`（本身按 N 个单元给出，一般有界谱浴，ED 已验证）；
- 实时：**QuAPI `Δt` 谱路径**——其块值即单元级全宽积分（`_gⱼₖ` 矩形 /
  `_gⱼⱼₖₖ` 三角，2026-11-05 再定性实证），Weyl 单元核 = `[1:N, 1:N]`
  直接截取。离散谱（Dirac-δ / DiscreteSpectrum）由 QuAPI 取**精确极点
  装配**（数值精确）；2026-11-07 起为唯一实时核路径（旧 Dirac-δ 闭式
  浴迹特例已删除，删除前互验 ≤4.1e-17）。
"""
function hybridization(bath, c::AbstractWeylContour; band::Int = 1)
    band == 1 || throw(ArgumentError("单元级核按带构造；多带用 influence_functional 复制"))
    blocks = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    if c isa WeylImagTime
        corr = Δτ(bath; N = c.Nτ, δτ = c.dτ)
        D = corr.data                             # ImagCorrelationFunction 内的 CorrelationMatrix，N 个单元
        blocks[(:τ, :τ)] = Matrix{ComplexF64}(D[1:c.Nτ, 1:c.Nτ])
    elseif c isa WeylMixedTime
        blocks = mixed_kernel_quapi(bath, c)
    else
        # 单模 Dirac-δ 浴统一走 QuAPI 谱路径（QuAPI 对离散谱取精确极点
        # 装配，δ 上与旧闭式浴迹路径互验 ≤4.1e-17，2026-11-07 移除特例）
        blocks = realtime_kernel_quapi(bath, c)
    end
    return HybridizationKernel(c, blocks, Float64(α_weight(bath)))
end

"""
    mixed_kernel_quapi(bath, c::WeylMixedTime) -> Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}

混合（Kadanoff–Baym L 形）轮廓的单元级核：QuAPI `Δm` 谱积分（任意有界谱），
九个分支块（(+,+)、(+,−)、(+,τ)、(−,+)、(−,−)、(−,τ)、(τ,+)、(τ,−)、(τ,τ)），
实时索引截取 [1:Nt]、虚时索引 [1:Nτ]——与实时核的单元截取同一语义。
"""
function mixed_kernel_quapi(bath, c::WeylMixedTime)
    qm = Δm(bath; Nτ = c.Nτ, t = c.Nt * c.dt, Nt = c.Nt, δτ = c.dτ)
    blocks = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
        Mq = GTEMPO.QuAPI.branch(qm; b1 = b1, b2 = b2)
        ri(c2) = c2 === :τ ? (1:c.Nτ) : (1:c.Nt)
        blocks[(b1, b2)] = Matrix{ComplexF64}(Mq[ri(b1), ri(b2)])
    end
    return blocks
end


"""
    realtime_kernel_quapi(bath, c::WeylRealTime) -> Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}

单元级实时核（全部浴型统一路径）：QuAPI `Δt` 的谱积分按单元索引
`[1:N, 1:N]` 截取。一般有界谱内部 quadgk；Dirac-δ / DiscreteSpectrum
离散谱由 QuAPI 取**精确极点装配**（数值上精确，无求积误差）。
（2026-11-07 起为唯一实时核路径：旧 Dirac-δ 闭式浴迹特例
`realtime_kernel_trace` 已删除，删除前互验 ≤4.1e-17。）
"""
function realtime_kernel_quapi(bath, c::WeylRealTime)
    qt = Δt(bath; N = c.Nt, t = c.Nt * c.dt)
    blocks = Dict{Tuple{Symbol,Symbol},Matrix{ComplexF64}}()
    for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-))
        Mq = GTEMPO.QuAPI.branch(qt; b1 = b1, b2 = b2)
        size(Mq, 1) == size(Mq, 2) == c.Nt + 1 ||
            throw(ArgumentError("QuAPI Δt 块尺寸 $(size(Mq)) 与单元数不符（应为 $(c.Nt+1)² 节点阵）"))
        blocks[(b1, b2)] = Matrix{ComplexF64}(Mq[1:c.Nt, 1:c.Nt])
    end
    return blocks
end

# ---------------------------------------------------------------------------
# GTEMPO hybriddynamics 的核视图：partialif_hybrid 通过
#   corr[conj(row ā), conj(col a)] 读取（"a-ā 约定 x[a,ā]=Δ"），
# 普通费米浴：a-ā → +Δ，ā-a → −Δ，同 conj → 0；跨支负号已内嵌在块值里，
# 由 hybriddynamics 的读数规则补回物理块 —— 与 QuAPI RealCorrelationFunction
# 的 getindex 语义完全一致，只是下标直接用单元 1..N。
# ---------------------------------------------------------------------------
struct KernelCorrView <: QuAPI.AbstractCorrelationFunction
    kernel::HybridizationKernel
    band::Int
end

function Base.getindex(v::KernelCorrView, i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.band == v.band && j.band == v.band) || return zero(ComplexF64)
    (i.conj != j.conj) || return zero(ComplexF64)
    Δ = v.kernel.blocks[(branch(i), branch(j))]
    # 混合轮廓的端点顶点（rt Nt+1）不进 Weyl 单元：越界读数 = 零耦合
    (1 <= i.j <= size(Δ, 1) && 1 <= j.j <= size(Δ, 2)) || return zero(ComplexF64)
    val = Δ[i.j, j.j]
    # 跨支读数符号约定（与 U1GTEMPO2 核值逐位核对，2026-11-05）：
    # ① 实时 (:+,:)：partialif_hybrid 消费端内嵌跨支负号（keldysh gate），
    #   视图预翻与之相消 ⟹ 净读数 = 核原始值；
    # ② 混合轮廓：消费端 gate 关闭；QuAPI Δm 的 τ↔rt 交叉块相对物理核差
    #   一个整体负号（ττ/±± 块无差——U1 对照实证）⟹ 视图仅翻 τ↔rt 读数；
    #   +↔− 读数沿用 ① 的双翻相消。
    if v.kernel.contour isa WeylRealTime
        # 实时：消费端 keldysh gate 内嵌跨支负号，视图预翻相消 ⟹ 净 = 原始值
        branch(i) !== branch(j) && (val = -val)
    elseif v.kernel.contour isa WeylMixedTime
        # 混合：仅翻 τ↔rt（Δm 交叉块相对物理核差整体负号——U1 核值对照实证）；
        # ±↔∓ 块与 ττ/±± 块直读原始值（符号组合扫描定案）。
        (branch(i) === :τ) !== (branch(j) === :τ) && (val = -val)
    end
    return i.conj ? -val : val
end

Z2Tensors.scalartype(::Type{KernelCorrView}) = ComplexF64

# KernelCorrView 的单元级块访问（扩展 QuAPI/GTEMPO 的 branch；勿新建同名函数，
# 否则会遮蔽 branch(::GrassmannContourIndex)）
import GTEMPO: branch
branch(x::KernelCorrView; b1::Symbol, b2::Symbol) = x.kernel.blocks[(b1, b2)]

