# influencefunctional/：IF 装配（plan/if_expmult_backend_plan.md W1 文件夹拆分
# + W3/W4 后端统一，2026-11-18；influence.jl 旧装配档案随入口切换退役删除）
#
#   partialif.jl   W3：fullfactor（Pairing/AnomalousPairing 整形式，紧凑站点表）+
#                  PartialIF（expmult!×PartialExp 幂零指数）+ PolarIF
#                  （expmult!×PolarExp SVD 秩一——同一算符第二种精确因子化）
#   exacttti.jl    W4 + Phase A′/B′：is_tti 门控 + 任意耦合模式（味对分组
#                  fit/direct）因子构造（ExpDecaySumSum 环族 / Pairing 容器）
#                  + ExactTTI（expmult!×PolarExp）；混合轴面退役
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
    PolarIF(; rtol=0.0, trunc=Defaults.alg_trunc_if())

PolarIF：fullfactor 的 SVD 秩一精确指数装配（expmult!×PolarExp）——同一
IF 算符 e^{Σ ā Δ a} 的第二种精确因子化（Pairing = 对角剥离 + 余阵 SVD；
AnomalousPairing = 味扇区折叠逐扇区 SVD）。适用范围与 PartialIF 一致
（任意轮廓 / kernel / 对称性；BCS 反常腿味在核键上，同一装配函数）。
`rtol` = PolarExp SVD 秩截断（相对 Frobenius 尾部² 口径；0 = 全数值秩，
机器精度精确——AnomalousPairing 跨扇区全局归一，BCS 配 rtol>0 时弱扇区
可能整扇区截掉）；`trunc` = 乘法压缩截断策略（Defaults.alg_trunc_if() =
truncdimcutoff 缺省）。
"""
@kwdef struct PolarIF{Tr<:TruncationStrategy} <: InfluenceFunctionalAlgorithm
    rtol::Float64 = 0.0
    trunc::Tr = Defaults.alg_trunc_if()
end

"""
    ExactTTI(; trunc=Defaults.alg_trunc_if(), alg_expan=Defaults.alg_expan())

ExactTTI{Tr<:TruncationStrategy,AE<:ExpExp.ExponentialExpansionAlgorithm}：
任意耦合模式核的 TTI 因子装配（2026-11-19 Phase A′/B′：is_tti 门控 +
band/味守卫退役；混合时间轴面退役——`tti_factors` 按轴 dispatch
（虚轴/实轴两方法，混合格无方法 = MethodError））。
`tti_factors` 按味对分组：虚时环族 ESS 极点压缩 / 实时 (g,l) 列拟合重构，
统一经 expmult!×PolarExp 装配（Pairing = SVD 秩一精确指数；ExpDecay 族委托
TTIExp 逐通道物化——要求 GM ≥ expdecayexp-same-orbital-fix a7b39bb 同轨道
触碰块修正）。**严格 TTI 门**（2026-11-19 审查裁决）：块非 Toeplitz、虚时
热配对 / 实时 Keldysh 审计失败一律 ArgumentError——ExactTTI = 严格 TTI
算法，非 TTI 核走 PartialIF/PolarIF；唯一非 fit 出口 = 虚时退化拟合节点
（|q| ≈ 0）的 Pairing 全矩阵兜底。

- `alg_expan`（类型参数 AE）= 函数指数分解算法（ExpExp.ExponentialExpansionAlgorithm
  面，Defaults.alg_expan() = AnchoredMatrixPencil 缺省，与铅笔路线同参）：
  fit 模式的核列分解经 `fit_exponentials(alg_expan, y)` 分派；分解精度
  随 alg 配置（精确单极点核机器精度；连续谱短核列受小样本极限）。
- `trunc` = 乘法压缩截断策略（SVDCompression(trunc=trunc)）。

装配路由：band 对角 normal 核 → 逐带独立装配返回逐带元组（单带核恒就地
单链）；含跨带/反常键（跨带耦合 / BCS 轨道）→ 恒单链（PartialIF 先例）。
"""
@kwdef struct ExactTTI{Tr<:TruncationStrategy,AE<:ExpExp.ExponentialExpansionAlgorithm} <: InfluenceFunctionalAlgorithm
    trunc::Tr = Defaults.alg_trunc_if()
    alg_expan::AE = Defaults.alg_expan()
end

"ExactTTI 的乘法压缩策略（trunc 随算法对象）"
_mult(alg::ExactTTI) = SVDCompression(trunc=alg.trunc)

"""
    hybriddynamics!(gmps, lattice, kernel, alg=PartialIF()) -> (GMPS | Tuple{GMPS...}, diagnostics)

原地 IF 装配（装配范围由 `kernel.bands` 决定，band 不重复输入；腿味在
`ContourBlockIndex.conj` 上，无调用端味旗标）。**PartialIF/PolarIF = 恒就地
装配 `gmps` 为单链**（`fullfactor` 单因子覆盖全部变量；原样返回，`gmps`
就地改写）。ExactTTI = band 对角 normal 核逐带独立装配：单带核就地单链；
多带核逐带从 `copy(gmps)` 起步返回逐带元组；**含跨带/反常键 → 恒单链**
（味对因子单序列覆盖全部变量，就地改写）。乘法压缩截断随算法对象携带
（`alg.trunc`）。第二输出 `diagnostics::Dict{Symbol,Any}` = 装配诊断：
`:truncerr`（乘法压缩总截断误差，Float64——第二输入 err 的累计终值）+
通道审计键（热配对/Keldysh 残差、逐对拟合展开与 fit/direct 模式）+
`:algorithm`/`:truncation`。ExactTTI 仅支持虚轴/实轴格（混合格无 ExactTTI
方法——MethodError，PartialIF/PolarIF 承接），且拟合层因子恒复 ⟹ 实链在装配入口复化（ExactTTI
返回态恒 ComplexF64，不再与 `gmps` 别名）。
"""
function hybriddynamics!(gmps, lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                         alg::InfluenceFunctionalAlgorithm=PartialIF())
    all(b -> 1 <= b <= lattice.bands, kernel.bands) ||
        throw(BoundsError(1:lattice.bands, kernel.bands))
    diagnostics = Dict{Symbol,Any}(:algorithm => alg, :truncation => alg.trunc)
    err = 0.0
    if alg isa PartialIF
        state, err = _partialif_state!(gmps, err, lattice, kernel, alg)
    elseif alg isa PolarIF
        state, err = _polarif_state!(gmps, err, lattice, kernel, alg)
    elseif alg isa ExactTTI
        state, err, diag = _exacttti_state!(gmps, err, lattice, kernel, alg)
        merge!(diagnostics, diag)
    else
        throw(ArgumentError("unsupported IF algorithm $(typeof(alg))"))
    end
    diagnostics[:truncerr] = err
    return state, diagnostics
end

"""
    hybriddynamics(lattice, kernel, alg=PartialIF()) -> (GMPS | Tuple{GMPS...}, diagnostics)

IF 链装配（真空起始）：`hybriddynamics(lattice, ...) =
hybriddynamics!(vacuumstate(promote_type(scalartype(lattice), scalartype(kernel)),
lattice), lattice, ...)`——提升的是两者的 scalartype。PartialIF/PolarIF 与
ExactTTI 跨带/反常核返回单 GMPS（恒单链）；ExactTTI band 对角 normal 多带核
返回逐带 `Tuple`（`Is...` 变参兼容——消费方对单 GMPS 须自包装 `(Is,)`）；
第二输出 = 装配 diagnostics（见 `hybriddynamics!`）。
"""
function hybriddynamics(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel,
                        alg::InfluenceFunctionalAlgorithm=PartialIF())
    return hybriddynamics!(vacuumstate(promote_type(TK.scalartype(lattice), TK.scalartype(kernel)), lattice),
                           lattice, kernel, alg)
end


include("partialif.jl")
include("exacttti.jl")
