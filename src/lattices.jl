# lattice 层：排布固定 Chrono（功能调整 #1；2026-11-16 起 ordering 类型层
# AbstractChronoOrdering 及其三个单例亦删除——排布是格点几何本身，不再是
# 可插拔策略，lattice 构造无 ordering 输入）。唯一排布的块序：
#   实时：[(+,1)…(+,Nt)][(−,Nt)…(−,1)]，K 时序链本征序直铺零 permute!；
#         每格 2M 位、(a,ā) 奇偶相邻 ⟹ AC 积分。
#   混合：[闲置 2 块][τ 1..Nτ][+ 1..Nt][− Nt…1]（闲置头块由真空穿透吸收）。
#   虚时：[格 Nτ][…][格 1]（格 i ⟷ 块 k−i+1）。
# 格点计数语义（与旧 lattice 的 k=步数+1 不同——Chrono 无额外端点格点）：
#   N/Nt/Nτ kwarg = 物理步数（= 单元数）；k/kt/kτ = 对应分支的格点数。
# lattice 按对称性单例类型参数化（2026-11-16：Symbol 字符串手选退役；单例
# 定义上游 GM/grassmannspaces.jl，本包再导出）——
# `symmetry::Sym` 字段携带 GM 的 GrassmannSymmetry 单例（Z2/Z2Z2/U1/U1U1），
# `TK.sectortype(lattice)` 编译期即得扇区类型（spec），下游免重复输入扇区
# 信息。ordering 管链排布、symmetry 管扇区标记，两者正交（§3.4）。

abstract type AbstractGrassmannLattice{Sym<:GrassmannSymmetry} end
abstract type AbstractImagGrassmannLattice{Sym} <: AbstractGrassmannLattice{Sym} end
abstract type AbstractRealGrassmannLattice{Sym} <: AbstractGrassmannLattice{Sym} end
abstract type AbstractMixedGrassmannLattice{Sym} <: AbstractGrassmannLattice{Sym} end

TK.sectortype(::Type{<:AbstractGrassmannLattice{Sym}}) where {Sym<:GrassmannSymmetry} = TK.sectortype(Sym)
TK.sectortype(lattice::AbstractGrassmannLattice) = TK.sectortype(lattice.symmetry)

branches(::Type{<:AbstractImagGrassmannLattice})  = (:τ,)
branches(::Type{<:AbstractRealGrassmannLattice})  = (:+, :-)
branches(::Type{<:AbstractMixedGrassmannLattice}) = (:τ, :+, :-)
branches(lattice::AbstractGrassmannLattice) = branches(typeof(lattice))

_index(x::AbstractGrassmannLattice, args...; kwargs...) =
    error("index not implemented for grassmann lattice type $(typeof(x))")
Base.getindex(x::AbstractGrassmannLattice, a::GrassmannContourIndex) =
    _index(x, a.j, conj=a.conj, branch=a.branch, band=a.band)

function indexmappings(lattice::AbstractGrassmannLattice)
    r = Dict{Tuple{Int, Bool, Symbol, Int}, Int}()
    for b in branches(lattice), i in 1:(b === :τ ? lattice.kτ : lattice.kt),
        c in (true, false), band in 1:lattice.bands
        r[(i, c, b, band)] = lattice[GrassmannContourIndex(i, b, band, c)]
    end
    return r
end

# ---------------------------------------------------------------------------
# 虚时
# ---------------------------------------------------------------------------

struct ImagGrassmannLattice{Sym} <: AbstractImagGrassmannLattice{Sym}
    δτ::Float64
    bands::Int
    N::Int              # τ 步数 = 格点数
    symmetry::Sym       # 对称性单例（Z2/Z2Z2/U1/U1U1）
    function ImagGrassmannLattice(δτ::Real, bands::Int, N::Int, symmetry::GrassmannSymmetry)
        bands_guard(TK.sectortype(symmetry), bands)   # 荷分辨形要求偶数带（轨道交错）
        N >= 1 || throw(ArgumentError("N ≥ 1 required"))
        new{typeof(symmetry)}(convert(Float64, δτ), bands, N, symmetry)
    end
end

function ImagGrassmannLattice(; δτ::Real, N::Int, bands::Int=1, symmetry::GrassmannSymmetry=Z2Symmetry())
    ImagGrassmannLattice(δτ, bands, N, symmetry)
end
Base.similar(x::ImagGrassmannLattice; δτ::Real=x.δτ, bands::Int=x.bands, N::Int=x.N,
             symmetry::GrassmannSymmetry=x.symmetry) =
    ImagGrassmannLattice(δτ, bands, N, symmetry)

function Base.getproperty(x::ImagGrassmannLattice, s::Symbol)
    if (s == :k) || (s == :kτ)
        return x.N
    elseif s == :τs
        return 0:x.δτ:x.N*x.δτ
    elseif s == :β
        return x.N * x.δτ
    elseif s == :T
        return 1 / x.β
    elseif s == :Nτ
        return x.N
    else
        getfield(x, s)
    end
end
TK.scalartype(::Type{<:AbstractImagGrassmannLattice}) = ComplexF64
Base.length(x::AbstractImagGrassmannLattice) = 2 * x.bands * x.k   # 每格 (a,ā)×bands

"虚时 Chrono _index：格 i（1..kτ）⟷ 链块 k−i+1（格 Nτ = 块 1）；
块内 (band 1: a,ā), …, (band M: a,ā)"
function _index(x::AbstractImagGrassmannLattice, i::Int; conj::Bool, band::Int=1, branch::Symbol=:τ)
    @boundscheck begin
        (branch === :τ) || throw(ArgumentError("branch must be :τ"))
        (1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
        (1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
    end
    block = x.k - i + 1
    return (block - 1) * 2x.bands + 2 * band - 1 + Int(conj)
end

# ---------------------------------------------------------------------------
# 实时（Keldysh）
# ---------------------------------------------------------------------------

struct RealGrassmannLattice{Sym} <: AbstractRealGrassmannLattice{Sym}
    δt::Float64
    bands::Int
    N::Int              # 实时步数 = 格点数（Chrono 无边界格点：ρ₀ 折进 (+,1)、闭合折进 (−,1)）
    symmetry::Sym       # 对称性单例
    function RealGrassmannLattice(δt::Real, bands::Int, N::Int, symmetry::GrassmannSymmetry)
        bands_guard(TK.sectortype(symmetry), bands)
        N >= 1 || throw(ArgumentError("N ≥ 1 required"))
        new{typeof(symmetry)}(convert(Float64, δt), bands, N, symmetry)
    end
end

RealGrassmannLattice(; δt::Real, N::Int, bands::Int=1, symmetry::GrassmannSymmetry=Z2Symmetry()) =
    RealGrassmannLattice(δt, bands, N, symmetry)
Base.similar(x::RealGrassmannLattice; δt::Real=x.δt, bands::Int=x.bands, N::Int=x.N,
             symmetry::GrassmannSymmetry=x.symmetry) =
    RealGrassmannLattice(δt, bands, N, symmetry)

function Base.getproperty(x::RealGrassmannLattice, s::Symbol)
    if (s == :k) || (s == :kt) || (s == :Nt)
        return x.N
    elseif s == :t
        return x.N * x.δt
    elseif s == :ts
        return 0:x.δt:x.N*x.δt
    else
        getfield(x, s)
    end
end
TK.scalartype(::Type{<:AbstractRealGrassmannLattice}) = ComplexF64
Base.length(x::AbstractRealGrassmannLattice) = 4 * x.bands * x.k   # 2k 块 × 2M 位

"实时 Chrono _index：块 = (+,i) → i；(−,i) → 2k−i+1；块内 2M 位"
function _index(x::AbstractRealGrassmannLattice, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
    @boundscheck begin
        (1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
        (1 <= i <= x.k) || throw(BoundsError(1:x.k, i))
        (branch in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    end
    block = branch === :+ ? i : 2x.k - i + 1
    return (block - 1) * 2x.bands + 2 * band - 1 + Int(conj)
end

# ---------------------------------------------------------------------------
# 混合（Kadanoff–Baym）
# ---------------------------------------------------------------------------

struct MixedGrassmannLattice{Sym} <: AbstractMixedGrassmannLattice{Sym}
    δt::Float64
    Nt::Int
    δτ::Float64
    Nτ::Int
    bands::Int
    symmetry::Sym       # 对称性单例
    function MixedGrassmannLattice(δt::Real, Nt::Int, δτ::Real, Nτ::Int, bands::Int, symmetry::GrassmannSymmetry)
        bands_guard(TK.sectortype(symmetry), bands)
        (Nt >= 1 && Nτ >= 1) || throw(ArgumentError("Nt, Nτ ≥ 1 required"))
        new{typeof(symmetry)}(convert(Float64, δt), Nt, convert(Float64, δτ), Nτ, bands, symmetry)
    end
end

function MixedGrassmannLattice(; δt::Real, Nt::Int, δτ::Real, Nτ::Int, bands::Int=1, symmetry::GrassmannSymmetry=Z2Symmetry())
    MixedGrassmannLattice(δt, Nt, δτ, Nτ, bands, symmetry)
end
Base.similar(x::MixedGrassmannLattice; δt::Real=x.δt, Nt::Int=x.Nt, δτ::Real=x.δτ, Nτ::Int=x.Nτ,
             bands::Int=x.bands, symmetry::GrassmannSymmetry=x.symmetry) =
    MixedGrassmannLattice(δt, Nt, δτ, Nτ, bands, symmetry)

function Base.getproperty(x::MixedGrassmannLattice, s::Symbol)
    if s == :t
        return x.Nt * x.δt
    elseif s == :β
        return x.Nτ * x.δτ
    elseif s == :T
        return 1 / x.β
    elseif s == :kt
        return x.Nt
    elseif s == :kτ
        return x.Nτ
    elseif s == :ts
        return 0:x.δt:x.t
    elseif s == :τs
        return 0:x.δτ:x.β
    else
        getfield(x, s)
    end
end
TK.scalartype(::Type{<:AbstractMixedGrassmannLattice}) = ComplexF64
Base.length(x::AbstractMixedGrassmannLattice) = 2 * x.bands * (x.Nτ + 2x.Nt)   # τ 格 + ± 格

"混合 Chrono _index（闲置头块已由真空穿透吸收、块号重排）：τ 格 j → 块 j；
+ 格 j → 块 Nτ+j；− 格 j → 块 Nτ+Nt+(Nt−j+1)"
function _index(x::AbstractMixedGrassmannLattice, i::Int; conj::Bool, branch::Symbol=:+, band::Int=1)
    @boundscheck begin
        (1 <= band <= x.bands) || throw(BoundsError(1:x.bands, band))
        if branch === :τ
            (1 <= i <= x.kτ) || throw(BoundsError(1:x.kτ, i))
        else
            (branch in (:+, :-)) || throw(ArgumentError("branch must be :+, :- or :τ"))
            (1 <= i <= x.kt) || throw(BoundsError(1:x.kt, i))
        end
    end
    block = branch === :τ ? i :
            branch === :+ ? x.Nτ + i :
            x.Nτ + x.Nt + (x.Nt - i + 1)
    return (block - 1) * 2x.bands + 2 * band - 1 + Int(conj)
end

# ---------------------------------------------------------------------------
# 工厂与构造约定（cell automaton / IF 构造的链几何视图）
# ---------------------------------------------------------------------------

"""
    GrassmannLattice(; contour::Symbol, kwargs...)

contour = :real/:Keldysh → RealGrassmannLattice；:imag → ImagGrassmannLattice；
:mixed/:Kadanoff → MixedGrassmannLattice。kwargs 透传（δt/Nt/δτ/Nτ/bands/
symmetry；无 ordering 输入——排布固定 Chrono，是格点几何本身）。
"""
function GrassmannLattice(; contour::Symbol, kwargs...)
    (contour in (:real, :imag, :Keldysh, :mixed, :Kadanoff)) ||
        throw(ArgumentError("contour must be :real (equivalent to :Keldysh), :imag or :mixed (equivalent to :Kadanoff)"))
    if (contour == :real) || (contour == :Keldysh)
        return RealGrassmannLattice(; kwargs...)
    elseif contour == :imag
        return ImagGrassmannLattice(; kwargs...)
    else
        return MixedGrassmannLattice(; kwargs...)
    end
end

"链块数（real 轮廓 2k、imag k+1、mixed Nτ+2Nt+3）"
function nblocks(lattice::AbstractGrassmannLattice)
    lattice isa AbstractRealGrassmannLattice && return 2 * lattice.k
    lattice isa AbstractImagGrassmannLattice && return lattice.k
    return lattice.Nτ + 2 * lattice.Nt
end

"""
    chronological_cells(lattice) -> Vector{Tuple{Symbol,Int}}

K 构造的时序格序：[τ 1..Nτ][+ 1..Nt][− Nt..1]（每格 = 该物理时间步的完整
单元；格 (b,j) 落在链块 `cell_block(lattice, (b,j))`）。
"""
function chronological_cells(lattice::AbstractGrassmannLattice)
    cells = Tuple{Symbol,Int}[]
    :τ in branches(lattice) && append!(cells, [(:τ, j) for j in 1:lattice.Nτ])
    :+ in branches(lattice) && append!(cells, [(:+, j) for j in 1:lattice.Nt])
    :+ in branches(lattice) && append!(cells, [(:-, j) for j in lattice.Nt:-1:1])
    return cells
end

"""
    cell_block(lattice, cell) -> Int

格 (b, j) 的链块号（_index 的 block 公式）。
"""
function cell_block(lattice::AbstractGrassmannLattice, cell::Tuple{Symbol,Int})
    b, j = cell
    lattice isa AbstractRealGrassmannLattice && return b === :+ ? j : 2lattice.k - j + 1
    lattice isa AbstractImagGrassmannLattice && return lattice.k - j + 1
    return b === :τ ? j :
           b === :+ ? lattice.Nτ + j :
           lattice.Nτ + lattice.Nt + (lattice.Nt - j + 1)
end

"块 → 变量位起点（每块 2·bands 位）"
block_offset(lattice::AbstractGrassmannLattice, block::Int) = (block - 1) * 2 * lattice.bands + 1

"""
    site_spaces(lattice) -> Vector{ElementarySpace}

链站空间序列（pair 布局 GMPS 每站 = 相邻变量对 (2s−1, 2s)；站 s 的 band =
mod1(s, bands)）——bridge 的 vacuumstate/pairize/引擎共用此视图。
"""
function site_spaces(lattice::AbstractGrassmannLattice)
    spec = TK.sectortype(lattice)
    nstations = length(lattice) ÷ 2
    return [station_space(spec, mod1(s, lattice.bands)) for s in 1:nstations]
end

"真空链（GM 真空：左边界真空、键恒等穿透）"
vacuumstate(lattice::AbstractGrassmannLattice) = GM.vacuumstate(ComplexF64, site_spaces(lattice))

"lattice 携带的对称性单例（Z2/Z2Z2/U1/U1U1；export 已声明）"
symmetry(l::AbstractGrassmannLattice) = l.symmetry
