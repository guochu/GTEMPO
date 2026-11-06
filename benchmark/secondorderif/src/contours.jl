# ============================================================================
# contours.jl — Weyl 轮廓：每条分支恰好 N 个全宽源单元（无额外浴传播端点）
#
# 对齐 U1GTEMPO2 的 RealTime/ImagTime（"N intervals have exactly N bath source
# cells"），但这里把单元直接映射到 GTEMPO GrassmannLattice 的顶点上，从而完整
# 复用 hybriddynamics / integrate / 观测量管线：
#
#   实时：RealGrassmannLattice(δt, bands, N=Nt−1) → lattice.k = Nt 个顶点，
#         顶点 j ⟷ 单元 j，每个顶点携带 (a_j⁺, ā_j⁺, a_j⁻, ā_j⁻)。
#   虚时：ImagGrassmannLattice(δτ, bands, N=Nτ) → lattice.kτ = Nτ+1 个顶点，
#         顶点 1 是 t=0 边界（平凡），单元 j 的 (a_j, ā_j) 放在顶点 j+1 ——
#         与 GTEMPO 一阶管线中 QuAPI 虚时核的 _latticeindex 位移约定一致。
# ============================================================================
abstract type AbstractWeylContour end

struct WeylRealTime <: AbstractWeylContour
    Nt::Int
    dt::Float64
    bands::Int
    function WeylRealTime(; Nt::Integer, dt::Real, bands::Int = 1)
        (Nt >= 1 && dt > 0 && isfinite(dt) && bands >= 1) ||
            throw(ArgumentError("Nt ≥ 1, dt > 0, bands ≥ 1 required"))
        new(Nt, Float64(dt), bands)
    end
end

struct WeylImagTime <: AbstractWeylContour
    Nτ::Int
    dτ::Float64
    bands::Int
    function WeylImagTime(; Nτ::Integer, dτ::Real, bands::Int = 1)
        (Nτ >= 1 && dτ > 0 && isfinite(dτ) && bands >= 1) ||
            throw(ArgumentError("Nτ ≥ 1, dτ > 0, bands ≥ 1 required"))
        new(Nτ, Float64(dτ), bands)
    end
end

nreal(c::WeylImagTime) = 0
nreal(c::WeylRealTime) = c.Nt
nimag(c::WeylRealTime) = 0
nimag(c::WeylImagTime) = c.Nτ
branches(::WeylRealTime) = (:+, :-)
branches(::WeylImagTime) = (:τ,)
ncells(c::AbstractWeylContour) = nreal(c) + nimag(c)

"单步复步长 λ：虚时 δτ、正向 +iδt、反向 −iδt（U1GTEMPO2 的 step(c,b)）"
stepλ(c::AbstractWeylContour, b::Symbol) =
    b === :τ ? c.dτ : b === :+ ? im * c.dt : -im * c.dt

"单元 j 覆盖的物理时间区间端点（t₀ = (j−1)·h, t₁ = j·h）"
cell_interval(c::WeylRealTime, j::Integer) = ((j - 1) * c.dt, j * c.dt)
cell_interval(c::WeylImagTime, j::Integer) = ((j - 1) * c.dτ, j * c.dτ)

"观察边界时刻 0..N（整数边界，U1GTEMPO2 的 times）"
boundary_times(c::WeylRealTime) = collect(0:c.Nt) .* c.dt
boundary_times(c::WeylImagTime) = collect(0:c.Nτ) .* c.dτ

"""
    ChronoKeldyshOrdering — 实时轮廓的时序顺序存储（完整版模块专用）

站点顺序 = K 时序链的本征顺序：[(+,1)…(+,N),(−,N)…(−,1)]，每块内
(band 1: a,ā), …, (band M: a,ā)。`_system_factor_rt` 的单元张量因此**原样
直铺**（时序 MPS 逐位不动；真空开边界）——旧路径经 permute!（swap 门链）
重排到顶点交错序，M≥2、N≥4 时中间键组合爆炸（N=4 >9 min）。

每格 2M 个位置、(a,ā) 从奇位开始 ⟹ AC 积分器的位置对 (2j−1,2j) 恰为共轭
对且永不跨格 ⟹ AdjacentConjugation：`_ac_integrate` 直接消费本格点。
"""
struct ChronoKeldyshOrdering <: GTEMPO.RealGrassmannOrdering end
GTEMPO.ConjugationStyle(::Type{ChronoKeldyshOrdering}) = GTEMPO.AdjacentConjugation()
GTEMPO.LayoutStyle(::Type{ChronoKeldyshOrdering}) = GTEMPO.TimeLocalLayout()
GTEMPO.RealTimeOrderingStyle(::Type{ChronoKeldyshOrdering}) = GTEMPO.TimeAscending()

# 顶点 i ∈ 1..k（k = 单元数），b = + → 块 i；b = − → 块 2k−i+1
function GTEMPO._index(x::AbstractRealGrassmannLattice{<:ChronoKeldyshOrdering}, i::Int;
                       conj::Bool, branch::Symbol = :+, band::Int = 1)
    @boundscheck begin
        (1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
        (1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
        (branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    end
    block = branch === :+ ? i : 2x.k - i + 1
    return (block - 1) * 2x.bands + 2 * band - 1 + Int(conj)
end

"""
    weyl_lattice(c::AbstractWeylContour; ordering=nothing) -> GrassmannLattice

把 Weyl 单元映射到 GTEMPO 格点（见文件头注释）。返回的格点满足：
实时 `lattice[k(ci)]` 给出单元变量在 MPS 中的 site 位置；虚时单元 j 位于
格点顶点 j+1。实时默认 `ChronoKeldyshOrdering`（时序顺序存储，K 免 permute!）。
"""
function weyl_lattice(c::WeylRealTime; ordering = nothing)
    ord = isnothing(ordering) ? ChronoKeldyshOrdering() : ordering
    # N = Nt−1 ⟹ k = Nt 个顶点 = Nt 个单元
    return RealGrassmannLattice(δt = c.dt, bands = c.bands, N = c.Nt - 1, ordering = ord)
end
function weyl_lattice(c::WeylImagTime; ordering = nothing)
    ord = isnothing(ordering) ? GTEMPO.A1Ā1B1B̄1() : ordering
    # N = Nτ ⟹ kτ = Nτ+1 个顶点；单元 j 在顶点 j+1（顶点 1 为 t=0 边界）
    return ImagGrassmannLattice(δτ = c.dτ, bands = c.bands, N = c.Nτ, ordering = ord)
end

"单元 j 的 (a, ā) 变量所在的格点顶点编号"
cell_vertex(c::WeylRealTime, j::Integer) = Int(j)
cell_vertex(c::WeylImagTime, j::Integer) = Int(j) + 1

"单元 j、分支 b 的 (a, ā) 一对 GrassmannContourIndex"
function cell_indices(c::AbstractWeylContour, j::Integer, b::Symbol, band::Int = 1)
    v = cell_vertex(c, j)
    return (GrassmannContourIndex(v, b, band, false),
            GrassmannContourIndex(v, b, band, true))
end
