# influencefunctional/：IF 装配（plan/if_expmult_backend_plan.md W1 文件夹拆分
# + W3/W4 后端统一，2026-11-18；influence.jl 旧装配档案随入口切换退役删除）
#
#   partialif.jl   W3：fullfactor（AnomalousPairing 整形式，紧凑站点表）+
#                  PartialIF（expmult!×PartialExp 幂零指数——语义换代，同一
#                  算符的第二种精确表示）
#   exacttti.jl    W4：tti_imag/real/cross_factors 三构造器（ExpDecaySumSum /
#                  Pairing 容器）+ ExactTTI（expmult!×PolarExp SVD 秩一）
#
# GM 前置 = polarexp-svd-rankone 计划（已实施）+ expdecayexp-same-orbital-fix
# （a7b39bb：ExpDecayExp 同轨道触碰块修正，ESS 通道严格性依赖）。装配范围由
# kernel 决定（band 信息在 `kernel.bands`，接口不再重复输入）：含跨带键的
# kernel（跨带耦合 / BCS 轨道）= 单链 IF，band 对角多带核 = 逐带 IF 元组。
# 乘法压缩截断随算法对象携带（PartialIF/ExactTTI.trunc = Defaults.alg_trunc_if()
# 缺省）；ExactTTI 的函数指数分解算法随 alg_expan 字段（Defaults.alg_expan()
# = AnchoredMatrixPencil 缺省）。

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
  Defaults.alg_expan() = AnchoredMatrixPencil 缺省，与铅笔路线同参）：tti_imag/
  real 的核列拟合经 `fit_exponentials(alg_expan, y)` 分派；分解精度随 alg 配置
  （精确单极点核机器精度；连续谱短核列受小样本极限）。
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
    hybriddynamics!(gmps, lattice, kernel, alg=PartialIF(); bcs=false)
        -> (GMPS | Tuple{GMPS...}, diagnostics)

原地 IF 装配（装配范围由 `kernel.bands` 决定，band 不重复输入）：含跨带键的
kernel（跨带耦合 / BCS 轨道）= 二次型不可按带分解 ⟹ 就地装配 `gmps` 为单链
（原样返回）；band 对角多带核 = 逐带从 `copy(gmps)` 起步装配（`gmps` 语义上
应为真空起始态，本身不变），返回逐带元组。乘法压缩截断随算法对象携带
（`alg.trunc`）。`bcs` = 跨带键的腿味解释（false = 一般法向有向耦合，true =
BCS 反常味表，见 `fullfactor`）。第二输出 `diagnostics::Dict{Symbol,Any}` =
装配诊断：`:truncerrs`（逐步乘法压缩的 GM truncerr）+ 通道审计键（Keldysh/
热配对残差、拟合展开）+ `:algorithm`/`:truncation`。ExactTTI 仅支持单带
normal 核（D7：BCS 仅 PartialIF）。
"""
function hybriddynamics!(gmps, lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                         alg::InfluenceFunctionalAlgorithm=PartialIF(); bcs::Bool=false)
    all(b -> 1 <= b <= lattice.bands, kernel.bands) ||
        throw(BoundsError(1:lattice.bands, kernel.bands))
    hascross = _has_cross_band_keys(kernel)
    diagnostics = Dict{Symbol,Any}(:algorithm => alg, :truncation => alg.trunc)
    if alg isa PartialIF
        state, truncerrs = _partialif_state(gmps, lattice, kernel, alg, bcs)
    elseif alg isa ExactTTI
        (bcs || hascross) && throw(ArgumentError(
            "ExactTTI does not support BCS/anomalous kernels (D7: BCS is PartialIF-only), " *
            "got bands $(kernel.bands)"))
        length(kernel.bands) == 1 || throw(ArgumentError(
            "ExactTTI requires a single-band normal kernel, got bands $(kernel.bands); " *
            "build per-band kernels for multi-band assembly"))
        state, truncerrs, diag = _exacttti_state(gmps, lattice, kernel, alg)
        merge!(diagnostics, diag)
    else
        throw(ArgumentError("unsupported IF algorithm $(typeof(alg))"))
    end
    diagnostics[:truncerrs] = truncerrs
    return state, diagnostics
end

"""
    hybriddynamics(lattice, kernel, alg=PartialIF(); bcs=false)
        -> (GMPS | Tuple{GMPS...}, diagnostics)

IF 链装配（真空起始）：`hybriddynamics(lattice, ...) =
hybriddynamics!(vacuumstate(promote_type(scalartype(lattice), scalartype(kernel)),
lattice), lattice, ...)`——提升的是两者的 scalartype。含跨带键的 kernel 返回
单 GMPS；band 对角多带 kernel 返回逐带 `Tuple`（`Is...` 变参兼容）；第二输出
= 装配 diagnostics（见 `hybriddynamics!`）。
"""
function hybriddynamics(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                        alg::InfluenceFunctionalAlgorithm=PartialIF(); bcs::Bool=false)
    return hybriddynamics!(vacuumstate(promote_type(TK.scalartype(lattice), TK.scalartype(kernel)), lattice),
                           lattice, kernel, alg; bcs)
end

"integrate(lattice, x0, x1...) —— Berezin 迹（链 = (K, Is...)，GM 引擎）"
integrate(lattice::AbstractGrassmannLattice, x0::GM.GMPS, x1::GM.GMPS...) = integrate(x0, x1...)
integrate_scaled(lattice::AbstractGrassmannLattice, x0::GM.GMPS, x1::GM.GMPS...) = integrate_scaled(x0, x1...)

include("partialif.jl")
include("exacttti.jl")
