# influencefunctional/：IF 装配（plan/if_expmult_backend_plan.md W1 文件夹拆分
# + W3/W4 后端统一，2026-11-18）
#
#   partialif.jl   W3：fullfactor（AnomalousPairing 整形式）+ PartialIF
#                  （expmult!×PartialExp 幂零指数——语义换代，同一算符的
#                  第二种精确表示，旧逐行秩一因子积冻结于 influence.jl）
#   exacttti.jl    W4：tti_imag/real/cross_factors 三构造器（ExpDecaySumSum /
#                  Pairing 容器）+ ExactTTI（expmult!×PolarExp SVD 秩一）
#
# GM 前置 = polarexp-svd-rankone 计划（PolarExp 支持 Pairing/AnomalousPairing
# SVD 秩一精确指数，已实施 GM main 789a1c6）。算法分派入口
# `influence_functional`（本文件）——旧 influence.jl 的 memory_cell/embed_cells
# 装配保留 deprecated 供每因子 dense 对拍（交叉验证期满退役，Phase 3/4）。

abstract type InfluenceFunctionalAlgorithm end

"""
    PartialIF(; maxdim=128, rtol=1e-12)

PartialIF：整形式 AnomalousPairing 幂零指数（PartialExp 截断 = 压缩）——
IF = exp(Σ ā Δ a) 的精确幂零级数表示（F6 语义换代：旧逐行秩一因子积的
同一算符）。SciPost Phys. Core 7, 063 (2024) 语义、GM 后端实现。
"""
Base.@kwdef struct PartialIF <: InfluenceFunctionalAlgorithm
    maxdim::Int = 128
    rtol::Float64 = 1e-12
end
_trunc(alg::PartialIF) = Truncation(D=alg.maxdim, rtol=alg.rtol)

"""
    ExactTTI(; rtol=1e-10, maxterms=40, cross_rtol=rtol, factor_rtol=nothing)

ExactTTI：tti_imag/real/cross 三通道因子（拟合照旧）统一经 expmult!×PolarExp
装配（Pairing/AnomalousPairing = SVD 秩一精确指数；ExpDecay 族委托 TTIExp
逐通道物化——要求 GM ≥ expdecayexp-same-orbital-fix a7b39bb 同轨道触碰块
修正）。real/imag 通道 SVD rtol = 0（精确）、cross 通道 = cross_rtol（截秩
口径同旧 cross 通道）。`factor_rtol` = 逐因子乘前压缩容差（nothing = 关闭）。
"""
struct ExactTTI <: InfluenceFunctionalAlgorithm
    rtol::Float64
    maxterms::Int
    cross_rtol::Float64
    factor_rtol::Union{Nothing,Float64}
    function ExactTTI(; rtol::Real=1e-10, maxterms::Integer=40, cross_rtol::Real=rtol,
                      factor_rtol::Union{Nothing,Real}=nothing)
        Float64(rtol); Int(maxterms)
        0 < rtol < 1 && 0 < cross_rtol < 1 && maxterms >= 1 ||
            throw(ArgumentError("invalid IF fit settings"))
        factor_rtol === nothing || (isfinite(factor_rtol) && 0 <= factor_rtol < 1) ||
            throw(ArgumentError("invalid factor cleanup tolerance"))
        new(Float64(rtol), Int(maxterms), Float64(cross_rtol),
            factor_rtol === nothing ? nothing : Float64(factor_rtol))
    end
end

struct InfluenceFunctional{K,S}
    kernel::K
    state::S
    band::Int
    diagnostics::Dict{Symbol,Any}
    reports::Vector{Float64}              # 逐步乘法压缩的 GM truncerr（标准标量通道）
end

"""
    influence_functional(lattice, kernel; alg=PartialIF(), band=1, trunc=Truncation())
        -> InfluenceFunctional

IF 构造入口（W3/W4 后端统一）：`band` = normal 带装配；BCS 轨道 kernel
（bands = [band1, band2]）走 fullfactor 轨道全域路径（`band` 取 kernel.bands
首元记账）。ExactTTI 不支持 BCS kernel（D7：BCS 仅 PartialIF）。
"""
function influence_functional(lattice::AbstractGrassmannLattice, kernel::HybridizationKernel;
                              alg::InfluenceFunctionalAlgorithm=PartialIF(), band::Int=1,
                              trunc::Truncation=Truncation())
    1 <= band <= lattice.bands || throw(BoundsError(1:lattice.bands, band))
    all(b -> 1 <= b <= lattice.bands, kernel.bands) || throw(BoundsError(1:lattice.bands, kernel.bands))
    isbcs = any(key[1].band != key[2].band for key in keys(kernel.blocks))
    band_eff = isbcs ? first(kernel.bands) : band
    isbcs || band in kernel.bands || throw(ArgumentError(
        "band $band not covered by kernel bands $(kernel.bands)"))
    diag = Dict{Symbol,Any}(:algorithm => alg, :truncation => trunc)
    state, truncerrs = if alg isa PartialIF
        _partialif_state(lattice, kernel; band=isbcs ? nothing : band, trunc)
    elseif alg isa ExactTTI
        isbcs && throw(ArgumentError(
            "ExactTTI does not support BCS kernels (D7: BCS is PartialIF-only), got orbital $(kernel.bands)"))
        _exacttti_state(lattice, kernel, alg, band, trunc)
    else
        throw(ArgumentError("unsupported IF algorithm $(typeof(alg))"))
    end
    diag[:assembly] = alg isa PartialIF ? :fullfactor_partial_exp : :tti_factors_polar_exp
    return InfluenceFunctional(kernel, state, band_eff, diag, truncerrs)
end

include("partialif.jl")
include("exacttti.jl")
