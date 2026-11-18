# influencefunctional/：IF 装配（plan/if_expmult_backend_plan.md W1 文件夹拆分
# + W3/W4 后端统一，2026-11-18；influence.jl 旧装配档案随入口切换退役删除）
#
#   partialif.jl   W3：fullfactor（AnomalousPairing 整形式）+ PartialIF
#                  （expmult!×PartialExp 幂零指数——语义换代，同一算符的
#                  第二种精确表示）
#   exacttti.jl    W4：tti_imag/real/cross_factors 三构造器（ExpDecaySumSum /
#                  Pairing 容器）+ ExactTTI（expmult!×PolarExp SVD 秩一）
#
# GM 前置 = polarexp-svd-rankone 计划（已实施）+ expdecayexp-same-orbital-fix
# （a7b39bb：ExpDecayExp 同轨道触碰块修正，ESS 通道严格性依赖）。装配范围由
# kernel 决定（band 信息在 `kernel.bands`，接口不再重复输入）：单带 normal /
# BCS 轨道核 = 单 IF，多带 normal 核 = 逐带 IF 元组。乘法压缩截断随算法对象
# 携带（PartialIF/ExactTTI.trunc = Defaults.alg_trunc_if() 缺省）；ExactTTI 的
# 函数指数分解算法随 alg_expan 字段（Defaults.alg_expan() 缺省）。

abstract type InfluenceFunctionalAlgorithm end

"""
    PartialIF(; trunc=Defaults.alg_trunc_if())

PartialIF：整形式 AnomalousPairing 幂零指数（PartialExp 截断 = 压缩）——
IF = exp(Σ ā Δ a) 的精确幂零级数表示（F6 语义换代：旧逐行秩一因子积的
同一算符）。SciPost Phys. Core 7, 063 (2024) 语义、GM 后端实现。
`trunc` = 乘法压缩截断策略（Defaults.alg_trunc_if() = truncdimcutoff 缺省）。
"""
@kwdef struct PartialIF{Tr<:TruncationStrategy} <: InfluenceFunctionalAlgorithm
    trunc::Tr = Defaults.alg_trunc_if()
end

"""
    ExactTTI(; trunc=Defaults.alg_trunc_if(), alg_expan=Defaults.alg_expan(),
             cross_rtol=1e-10, factor_rtol=nothing)

ExactTTI：tti_imag/real/cross 三通道因子统一经 expmult!×PolarExp 装配
（Pairing/AnomalousPairing = SVD 秩一精确指数；ExpDecay 族委托 TTIExp 逐通道
物化——要求 GM ≥ expdecayexp-same-orbital-fix a7b39bb 同轨道触碰块修正）。

- `alg_expan` = 函数指数分解算法（ExpExp.ExponentialExpansionAlgorithm 面，
  Defaults.alg_expan() = OverDeterminedProny 缺省）：tti_imag/real 的核列
  拟合层经 `fit_exponentials(alg_expan, y)` 分派；拟合精度随算法配置
  （缺省 Prony(n=15, tol=1e-4) 对连续谱核列残差 ~1e-5–1e-7 档，精确单极点核
  机器精度；更高精度配 `Defaults.alg_expan(tol=1e-10, n=40)` 或 MatrixPencil）。
- `trunc` = 乘法压缩截断策略（SVDCompression(trunc=trunc)）；
- `cross_rtol` = cross 通道 Pairing SVD 截秩（real/imag 通道恒 rtol=0 精确）；
- `factor_rtol` = 逐因子乘前压缩容差（nothing = 用 trunc）。
"""
struct ExactTTI{Tr<:TruncationStrategy} <: InfluenceFunctionalAlgorithm
    trunc::Tr
    alg_expan::ExpExp.ExponentialExpansionAlgorithm
    cross_rtol::Float64
    factor_rtol::Union{Nothing,Float64}
    function ExactTTI(; trunc=Defaults.alg_trunc_if(), alg_expan=Defaults.alg_expan(),
                      cross_rtol::Real=1e-10, factor_rtol::Union{Nothing,Real}=nothing)
        isfinite(cross_rtol) && 0 < cross_rtol < 1 ||
            throw(ArgumentError("invalid cross_rtol $cross_rtol"))
        factor_rtol === nothing || (isfinite(factor_rtol) && 0 <= factor_rtol < 1) ||
            throw(ArgumentError("invalid factor cleanup tolerance"))
        new{typeof(trunc)}(trunc, alg_expan, Float64(cross_rtol),
            factor_rtol === nothing ? nothing : Float64(factor_rtol))
    end
end

"ExactTTI 的乘法压缩策略（trunc 随算法对象；factor_rtol 覆盖为纯 rtol 截断）"
_mult(alg::ExactTTI) =
    alg.factor_rtol === nothing ? SVDCompression(trunc=alg.trunc) :
    SVDCompression(trunc=Truncation(D=typemax(Int), rtol=alg.factor_rtol))

"GMPS 起始态深拷贝（多带核逐带装配的起步态；TensorMap 逐张拷贝）"
_copy_mps(x::GM.GMPS) = begin
    out = GM.GMPS(copy.(x.data))
    GM.setlogscale!(out, GM.logscale(x))
    out
end

"""
    hybriddynamics!(gmps, lattice, kernel, alg=PartialIF(); diagnostics=Dict()) -> GMPS | Tuple{GMPS...}

原地 IF 装配（装配范围由 `kernel.bands` 决定，band 不重复输入）：单带
normal / BCS 轨道核 = 就地装配 `gmps`（原样返回）；多带 normal 核 = 逐带从
`copy(gmps)` 起步装配（`gmps` 语义上应为真空起始态，本身不变），返回逐带
元组。乘法压缩截断随算法对象携带（`alg.trunc`）；truncerrs 与通道审计
包装进 `diagnostics`（`:truncerrs` + 审计键 + `:algorithm`/`:truncation`）。
ExactTTI 不支持 BCS kernel（D7：BCS 仅 PartialIF）。
"""
function hybriddynamics!(gmps, lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                         alg::InfluenceFunctionalAlgorithm=PartialIF();
                         diagnostics::Dict{Symbol,Any}=Dict{Symbol,Any}())
    all(b -> 1 <= b <= lattice.bands, kernel.bands) ||
        throw(BoundsError(1:lattice.bands, kernel.bands))
    isbcs = _is_orbital_kernel(kernel)
    diagnostics[:algorithm] = alg
    diagnostics[:truncation] = alg.trunc
    if alg isa PartialIF
        state, truncerrs = _partialif_state(gmps, lattice, kernel, alg)
    elseif alg isa ExactTTI
        isbcs && throw(ArgumentError(
            "ExactTTI does not support BCS kernels (D7: BCS is PartialIF-only), got orbital $(kernel.bands)"))
        state, truncerrs, diag = _exacttti_state(gmps, lattice, kernel, alg)
        merge!(diagnostics, diag)
    else
        throw(ArgumentError("unsupported IF algorithm $(typeof(alg))"))
    end
    diagnostics[:truncerrs] = truncerrs
    return state, diagnostics
end

"""
    hybriddynamics(lattice, kernel, alg=PartialIF(); diagnostics=Dict()) -> GMPS | Tuple{GMPS...}

IF 链装配（真空起始）：`hybriddynamics(lattice, ...) =
hybriddynamics!(vacuumstate(promote_type(TK.scalartype(lattice), TK.scalartype(kernel)),
lattice), lattice, ...)`——提升的是两者的 scalartype。单带 kernel 返回单
GMPS；多带 kernel 返回逐带 `Tuple`（`Is...` 变参兼容）。
"""
function hybriddynamics(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                        alg::InfluenceFunctionalAlgorithm=PartialIF();
                        diagnostics::Dict{Symbol,Any}=Dict{Symbol,Any}())
    return hybriddynamics!(vacuumstate(promote_type(TK.scalartype(lattice), TK.scalartype(kernel)), lattice),
                           lattice, kernel, alg; diagnostics)
end

"integrate(lattice, x0, x1...) —— Berezin 迹（链 = (K, Is...)，GM 引擎）"
integrate(lattice::AbstractGrassmannLattice, x0::GM.GMPS, x1::GM.GMPS...) = integrate(x0, x1...)
integrate_scaled(lattice::AbstractGrassmannLattice, x0::GM.GMPS, x1::GM.GMPS...) = integrate_scaled(x0, x1...)

include("partialif.jl")
include("exacttti.jl")
