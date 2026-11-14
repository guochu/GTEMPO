# bridge.jl：GrassmannMPS 后端适配层（对称性泛型化）
#
# 布局映射（plan/gtempo_grassmannmps_backend_plan.md §2，D1 裁决）：
#   构造路径（sysdynamics/influencefunctional）保持变量粒度（每变量一站、2 维
#   物理腿、bosonic 扇区），经 `pairize` 相邻同味站对楔积融合成 pair 布局 GMPS
#   （每胞 2 站：up 形/down 形 grassmannspace，4 维站）后进入 GrassmannMPS。
#
# 对称性用户手选（2026-11-13 D1 裁决，§3.4）：四个 spec 单例 ↔ :z2/:z2z2/:u1/:u1u1。
# 荷标签（构造层纯数据）类型随 spec：
#   :z2   Int（各带宇并合并——旧 GTEMPO 语义；sector = mod 2）
#   :z2z2 NTuple{2,Int}（每带宇并，bands≤2；sector = mod 2 积）
#   :u1   Int（各带总粒子数合并；sector = 原值）
#   :u1u1 NTuple{2,Int}（每带荷，bands≤2；U1GTEMPO2 同款）
# 费米符号不在标签里——由 GM.GrassmannBackend 显式注入（@grassmann 宏）。
# 参考实现 = U1GTEMPO2 bridge.jl（ChargeCharge 通道机器精度验证过）；
# 其余三通道靠 GM abelian 电池 + 本包跨对称性一致电池 + ED 三层兜底（R9）。

import GrassmannMPS as GM
import TensorKit as TK
using TensorKit: ⊗, ⊠, ←, SectorVector, SectorDict
import MatrixAlgebraKit as MAK
using LinearAlgebra: Diagonal, I, norm, diag

# ---------------------------------------------------------------------------
# 对称性 spec（扇区类型注册面：本文件是唯一允许按对称性分支的位置之一，
# 与 GM/grassmannspaces.jl 的注册表对表）
# ---------------------------------------------------------------------------

abstract type AbstractSymmetrySpec end
struct Z2Spec <: AbstractSymmetrySpec end       # :z2   全带合并宇并（旧语义）
struct Z2Z2Spec <: AbstractSymmetrySpec end     # :z2z2 band1/2 各占宇并槽（≤2 带）
struct U1Spec <: AbstractSymmetrySpec end       # :u1   全带合并总荷
struct U1U1Spec <: AbstractSymmetrySpec end     # :u1u1 band1→槽1 band2→槽2（≤2 带）

const SYMMETRY_SPECS = IdDict{Symbol,AbstractSymmetrySpec}(
    :z2 => Z2Spec(), :z2z2 => Z2Z2Spec(), :u1 => U1Spec(), :u1u1 => U1U1Spec())

"Symbol → spec（未知 Symbol 显式拒绝，禁静默降级——D1 裁决）"
symmetry_spec(s::Symbol) = get(SYMMETRY_SPECS, s, nothing) === nothing ?
    throw(ArgumentError("unknown symmetry :$s (allowed: :z2, :z2z2, :u1, :u1u1)")) :
    SYMMETRY_SPECS[s]

# 荷标签代数。
# U1 形（:u1/:u1u1）：整数荷（Int / NTuple{2,Int}），普通加减。
# 宇并形（:z2/:z2z2）：标签 = 奇并向量（Z2Label 包装 Int，算术 mod 2——
#   +1 与 −1 同扇区，±1 整数算术会把两扇区压重）；宇并加法 = 异或、
#   负号 = 恒等。普通泛型代码对标签只做 qadd/qsub/qneg/==/sort/hash。
struct Z2Label
    n::Int
    Z2Label(n::Integer) = new(mod(Int(n), 2))
end
Base.:(==)(a::Z2Label, b::Z2Label) = a.n == b.n
Base.hash(a::Z2Label, h::UInt) = hash(a.n, h)
Base.isless(a::Z2Label, b::Z2Label) = a.n < b.n
Base.show(io::IO, a::Z2Label) = print(io, "ζ", a.n)

@inline qadd(a::Int, b::Int) = a + b
@inline qsub(a::Int, b::Int) = a - b
@inline qneg(a::Int) = -a
@inline qadd(a::NTuple{2,Int}, b::NTuple{2,Int}) = (a[1] + b[1], a[2] + b[2])
@inline qsub(a::NTuple{2,Int}, b::NTuple{2,Int}) = (a[1] - b[1], a[2] - b[2])
@inline qneg(a::NTuple{2,Int}) = (-a[1], -a[2])
@inline qadd(a::Z2Label, b::Z2Label) = Z2Label(a.n + b.n)
@inline qsub(a::Z2Label, b::Z2Label) = Z2Label(a.n + b.n)
@inline qneg(a::Z2Label) = a
@inline qadd(a::NTuple{2,Z2Label}, b::NTuple{2,Z2Label}) = (qadd(a[1], b[1]), qadd(a[2], b[2]))
@inline qsub(a::NTuple{2,Z2Label}, b::NTuple{2,Z2Label}) = (qadd(a[1], b[1]), qadd(a[2], b[2]))
@inline qneg(a::NTuple{2,Z2Label}) = a

@inline zero_charge(::Z2Spec) = Z2Label(0)
@inline zero_charge(::U1Spec) = 0
@inline zero_charge(::Z2Z2Spec) = (Z2Label(0), Z2Label(0))
@inline zero_charge(::U1U1Spec) = (0, 0)

"荷标签 → bosonic 扇区（费米符号在 backend 符号层）"
sector(::Z2Spec, q::Z2Label) = TK.Z2Irrep(q.n)
sector(::U1Spec, q::Int) = TK.U1Irrep(q)
sector(::Z2Z2Spec, q::NTuple{2,Z2Label}) = Z2Z2((TK.Z2Irrep(q[1].n), TK.Z2Irrep(q[2].n)))
sector(::U1U1Spec, q::NTuple{2,Int}) = ChargeCharge(TK.U1Irrep(q[1]), TK.U1Irrep(q[2]))

"bosonic 扇区 → 荷标签（LegacyTruncation 排序与 _convert_site 用）"
charge_label(::Z2Spec, s::TK.Z2Irrep) = Z2Label(Int(s.n))
charge_label(::U1Spec, s::TK.U1Irrep) = Int(s.charge)
charge_label(::Z2Z2Spec, s::Z2Z2) = (Z2Label(Int(s[1].n)), Z2Label(Int(s[2].n)))
charge_label(::U1U1Spec, s::ChargeCharge) = (Int(s[1].charge), Int(s[2].charge))

"扇区类型（qspace 构造用）"
sector_type(::Z2Spec) = TK.Z2Irrep
sector_type(::Z2Z2Spec) = Z2Z2
sector_type(::U1Spec) = TK.U1Irrep
sector_type(::U1U1Spec) = ChargeCharge

"扇区类型 → spec 反查（链的键扇区类型即定 spec；链上操作无需显式传 spec）"
sector_spec(::Type{TK.Z2Irrep}) = Z2Spec()
sector_spec(::Type{TK.U1Irrep}) = U1Spec()
sector_spec(::Type{ChargeCharge}) = U1U1Spec()
sector_spec(::Type{Z2Z2}) = Z2Z2Spec()
chain_spec(x::GM.GMPS) = sector_spec(TK.sectortype(x[1]))

"荷分辨形 bands 守卫（:z2z2/:u1u1 每带分辨要求 bands≤2，超限显式拒绝——D1）"
function bands_guard(s::Symbol, bands::Int)
    if (s === :z2z2 || s === :u1u1) && bands > 2
        throw(ArgumentError(
            "symmetry :$s supports at most 2 bands (got $bands); use :z2 or :u1 for merged marks"))
    end
    return bands
end

"变量荷：band b 上 η 生成元的荷标签（η̄ = qneg）
——宇并形算术 mod 2 后 η̄ ≡ η（同一标签），带内 (a,ā) 区分在站形槽位，
不在荷上；荷分辨形 band→槽位映射照旧。"
var_charge(::Z2Spec, band::Int, bands::Int) = Z2Label(1)
var_charge(::U1Spec, band::Int, bands::Int) = 1
var_charge(::Z2Z2Spec, band::Int, bands::Int) = band == 1 ? (Z2Label(1), Z2Label(0)) : (Z2Label(0), Z2Label(1))
var_charge(::U1U1Spec, band::Int, bands::Int) = band == 1 ? (1, 0) : (0, 1)

"站形：band b 的物理站空间（单荷形全站同形；双荷形按带奇偶交替朝向——
:u1u1 = U1GTEMPO2 同款 up/down 交替站形）"
station_space(::Z2Spec, band::Int) = GM.grassmannspace(TK.Z2Irrep)
station_space(::U1Spec, band::Int) = GM.grassmannspace(TK.U1Irrep)
station_space(::Z2Z2Spec, band::Int) = GM.grassmannspace(Z2Z2, isodd(band) ? :up : :down)
station_space(::U1U1Spec, band::Int) = GM.grassmannspace(ChargeCharge, isodd(band) ? :up : :down)

"变量荷 → 站朝向（_convert_site 用；单荷形无朝向 → 返回 nothing 由调用方取全站同形）"
station_flavor(::Union{Z2Spec,U1Spec}, q) = :up   # 单荷形 spin 被忽略，恒 :up 合法
station_flavor(::Z2Z2Spec, q) = iszero(q[2].n) ? :up : :down   # Z2Label 须拆 .n 比较（== Int 恒 false）
station_flavor(::U1U1Spec, q) = q[2] == 0 ? :up : :down

"键标签序列 → 分级空间（标签计数聚合；排序 = 标签升序保证扇区插入序确定）"
function qspace(spec::AbstractSymmetrySpec, labels::AbstractVector)
    G = sector_type(spec)
    TK.Vect[G](sector(spec, q) => count(==(q), labels) for q in sort(unique(labels)))
end
"2 维物理腿空间：vacuum ⊕ generator q"
pspace(spec::AbstractSymmetrySpec, q) = qspace(spec, [zero_charge(spec), q])
"物理腿的非平凡扇区（2 维腿只有一个生成元扇区）"
occupied(p) = only(filter(!=(one(TK.sectortype(p))), collect(TK.sectors(p))))

"每胞变量荷序列（Chrono 排布每格 = [(band b: a, ā)] 逐带）"
function cell_charges(spec::AbstractSymmetrySpec, bands::Int)
    out = Any[]
    for b in 1:bands
        qc = var_charge(spec, b, bands)
        push!(out, qc, qneg(qc))
    end
    return out
end

"BCS 守卫（D7 裁决）：BCS 反常关联破坏粒子数守恒、宇称恒守恒 → 仅 :z2"
bcs_symmetry_guard(sym::Symbol) =
    sym === :z2 || throw(ArgumentError(
        "BCS bath supports only symmetry=:z2 (anomalous pairing breaks particle-number " *
        "conservation; only parity is conserved), got :$sym"))

# ---------------------------------------------------------------------------
# 链几何：变量 ↔ 站映射（Chrono 排布每格 2·bands 变量 = bands 站）
# ---------------------------------------------------------------------------

"链站空间序列：站 s 的 band = mod1(s, bands)"
function site_spaces(spec::AbstractSymmetrySpec, nstations::Int, bands::Int)
    [station_space(spec, mod1(s, bands)) for s in 1:nstations]
end

"变量位（1..2·nstations）→ 站位"
@inline station_of(varpos::Int) = (varpos + 1) ÷ 2

# ---------------------------------------------------------------------------
# 截断（U1GTEMPO2 D1/D2 裁决移植：有状态 TruncationStrategy 逐位复刻旧
# select_dimensions + CompressionReport 记账；荷排序键 = 标签自身，跨荷全局预算）
# ---------------------------------------------------------------------------

struct Truncation
    maxdim::Int
    rtol::Float64
    cutoff::Float64
    function Truncation(; maxdim::Integer=128, rtol::Real=1e-12, cutoff::Real=0)
        maxdim > 0 && isfinite(rtol) && 0 <= rtol < 1 && isfinite(cutoff) && cutoff >= 0 ||
            throw(ArgumentError("invalid truncation"))
        new(maxdim, rtol, cutoff)
    end
end

struct CompressionReport
    discarded_norms::Vector{Float64}
    relative_errors::Vector{Float64}
    bond_dimensions::Vector{Int}
    saturated::Bool
    discarded_logabs::Vector{Float64}
end
CompressionReport(d, r, b, s) = CompressionReport(d, r, b, s, log.(d))

"""
compress! 内部的有状态截断策略：`svd_trunc` 逐键调用 `findtruncated_svd`，
选择语义逐位复刻旧 `select_dimensions`（跨荷全局预算 + 相对尾容差 + 绝对
cutoff），并在调用序（键 L−1..1）里记账。荷排序键 = 荷标签自身
（Int / NTuple{2,Int} 皆可比较，跨荷预算语义不变）。
"""
mutable struct LegacyTruncation{L} <: MAK.TruncationStrategy
    spec::AbstractSymmetrySpec
    maxdim::Int
    rtol::Float64
    cutoff::Float64
    l0::Float64                       # 扫掠起点 logscale（QR 归一之后）
    dims::Vector{Dict{L,Int}}
    discarded::Vector{Float64}
    relative::Vector{Float64}
    saturated::Vector{Bool}
    spectral_norms::Vector{Float64}   # 每键 sqrt(Σ|s|²)（rightorth! 吸收进 logscale 的量）
end
LegacyTruncation(spec::AbstractSymmetrySpec, ::Type{L}, trunc::Truncation, l0::Float64) where L =
    LegacyTruncation{L}(spec, trunc.maxdim, trunc.rtol, iszero(trunc.cutoff) ? 0.0 : trunc.cutoff, l0,
                        Dict{L,Int}[], Float64[], Float64[], Bool[], Float64[])

function MAK.findtruncated_svd(values::TK.SectorVector, s::LegacyTruncation{L}) where L
    entries = sort!([(v, charge_label(s.spec, q), j) for (q, vs) in pairs(values) for (j, v) in enumerate(vs)];
                    by = x -> (-x[1], x[2], x[3]))
    tails = zeros(length(entries) + 1)
    for i in length(entries):-1:1
        tails[i] = tails[i+1] + abs2(entries[i][1])
    end
    nrm = sqrt(tails[1])
    # 谱归一漂移逐键累积：本次调用时链相对扫掠起点的 logscale 增量 = Σ_{j<本次} log nrm_j
    drift = sum(log(nz) for nz in s.spectral_norms if !iszero(nz); init=0.0)
    cutoff_eff = iszero(s.cutoff) ? 0.0 : s.cutoff * exp(-drift)
    k = 0
    for (v, _, _) in entries
        k >= s.maxdim && break
        k > 0 && (tails[k+1] <= s.rtol^2 * tails[1] || v <= cutoff_eff) && break
        k += 1
    end
    dims = Dict{L,Int}()
    for (_, q, j) in entries[1:k]
        dims[q] = max(get(dims, q, 0), j)
    end
    push!(s.dims, dims); push!(s.discarded, sqrt(tails[k+1]))
    push!(s.relative, iszero(tails[1]) ? 0.0 : sqrt(tails[k+1] / tails[1]))
    push!(s.saturated, k >= s.maxdim)
    push!(s.spectral_norms, nrm)
    return TK.SectorDict(q => (1:get(dims, charge_label(s.spec, q), 0)) for q in keys(values))
end

"结构完整性检查（旧 check_charge 语义）：非空、logscale 有限、左边界真空、
键匹配、物理站为规范同站形"
function _check_structure(x::GM.GMPS)
    isempty(x) && throw(ArgumentError("empty GMPS"))
    isfinite(GM.logscale(x)) || throw(ArgumentError("nonfinite GMPS logscale"))
    TK.space(x[1], 1) == oneunit(TK.space(x[1], 1)) || throw(ArgumentError("left boundary must be vacuum"))
    for i in 1:length(x)-1
        TK.domain(x[i])[1] == GM.space(x[i+1], 1) || throw(ArgumentError("bond mismatch between sites $i and $(i+1)"))
    end
    foreach(i -> GM.grassmann_basis(GM.physpace(x, i)), 1:length(x))
    true
end

"""
    compress!(x::GMPS; trunc=Truncation()) -> CompressionReport

QR 前扫（GM.leftorth!）+ SVD 回扫（GM.rightorth! 携带 LegacyTruncation）；
cutoff 相对原系数链尺度换算（logscale 漂移），报告的 discarded_logabs = 键
丢弃范数 + 扫掠后 logscale（旧实现逐位语义）。零链早退保持。
"""
function compress!(x::GM.GMPS; trunc::Truncation=Truncation())
    _check_structure(x)
    original_logscale = GM.logscale(x)
    # QR 前扫带逐步 rescale（旧 compress! 语义：乘积链站点幅值可跨越 1e100，
    # 必须在每步吸收防溢出级联——GM.leftorth! 只做入口一次性 rescale，不够）
    GM._rescale_sites!(x)
    for i in 1:length(x)-1
        n = norm(x[i])
        if !iszero(n) && (n > 1e100 || n < 1e-100)
            x.data[i] = x[i] / n
            GM.setlogscale!(x, GM.logscale(x) + log(n))
        end
        GM._left_step!(x, i)
    end
    n = norm(x[end])
    if !iszero(n) && (n > 1e100 || n < 1e-100)
        x.data[end] = x[end] / n
        GM.setlogscale!(x, GM.logscale(x) + log(n))
    end
    if iszero(norm(x[end]))
        for i in eachindex(x.data)
            x.data[i] = TK.zeros(TK.scalartype(x[i]), TK.codomain(x[i]) ← TK.domain(x[i]))
        end
        return CompressionReport(zeros(length(x) - 1), zeros(length(x) - 1), bond_dimensions(x), false)
    end
    # 键谱 = 满尺度值（rightorth! 在 svd_trunc 之后才吸收 Σs²）；cutoff 相对
    # 原链尺度 = 旧 cutoff_log 逻辑（漂移 = QR 扫掠后的 logscale 增量）
    cutoff_log = iszero(trunc.cutoff) ? -Inf : log(trunc.cutoff) - (GM.logscale(x) - original_logscale)
    cutoff = cutoff_log == -Inf ? 0.0 :
             cutoff_log > log(floatmax(Float64)) ? floatmax(Float64) : exp(cutoff_log)
    strategy = LegacyTruncation(chain_spec(x), typeof(zero_charge(chain_spec(x))),
                                Truncation(maxdim=trunc.maxdim, rtol=trunc.rtol, cutoff=cutoff), GM.logscale(x))
    GM.rightorth!(x; trunc=strategy)
    L = length(x)
    dn = zeros(L - 1); re = similar(dn)
    for (m, record) in enumerate(strategy.discarded)
        dn[L-m] = record; re[L-m] = strategy.relative[m]
    end
    # Δpoly（多项式尺度）= dn_m·exp(截断时刻 logscale)：bond m 丢弃尾巴的时刻，
    # 多项式 = exp(L_m)·C(tensors)，故换算 = log dn_m + L_m（逐键漂移累积；
    # L₀ 已在构造时传入 strategy——rightorth! 末站归一发生在丢弃之后，不进 L_m）
    run = strategy.l0
    logs = similar(dn)
    for m in 1:L-1
        bond = L - m
        logs[bond] = log(dn[bond]) + run
        iszero(strategy.spectral_norms[m]) || (run += log(strategy.spectral_norms[m]))
    end
    absolute = [v == -Inf ? 0.0 : v > log(floatmax(Float64)) ? Inf : exp(v) for v in logs]
    CompressionReport(absolute, re, bond_dimensions(x), any(strategy.saturated), logs)
end

# ---------------------------------------------------------------------------
# 乘法（GrassmannMPS 朴素外积 + compress!；old multiply 逐位语义）
# ---------------------------------------------------------------------------

function multiply(x::GM.GMPS, y::GM.GMPS; trunc::Union{Nothing,Truncation}=Truncation(), reports=nothing)
    out = x * y
    if trunc !== nothing
        r = compress!(out; trunc)
        reports === nothing || push!(reports, r)
    end
    out
end

# ---------------------------------------------------------------------------
# pairize：变量粒度裸张量链 → pair 布局 GMPS（对称性泛型版）
# ---------------------------------------------------------------------------

"站形（flavor）：物理生成元荷 → 站朝向（单荷形恒 :up）"
_station_space(spec::AbstractSymmetrySpec, q) = station_space(spec, station_flavor(spec, q) === :up ? 1 : 2)

"变量粒度站张量（cod=(键出,物理 2 维), dom=(键入)）→ pair 布局 4 维站张量：
**键腿恒等**（荷标签扇区原样保留——荷分辨对称性的效率即在于此，键两端恒等
映射故必一致）；物理腿按味嵌入各自站形的槽——vac→unit 槽 1、η 变量 →
grassmann_basis.eta 标签、η̄ 变量 → grassmann_basis.etabar 标签（异扇区形 =
扇区分辨、自对偶形 = 同扇区分槽——U1GTEMPO2 原版只支持异扇区，自对偶槽位
是本泛化版的关键补丁），pair 槽留空（逐对 _fuse_adjacent 楔积时填充）。
`conj` = 变量的共轭旗（pairize 内奇位 = a/η、偶位 = ā/η̄）。"
function _convert_site(spec::AbstractSymmetrySpec, t, conj::Bool)
    q = charge_label(spec, occupied(GM.space(t, 2)))
    E = _station_space(spec, q)
    basis = GM.grassmann_basis(E)
    tag = conj ? basis.etabar : basis.eta
    G = TK.sectortype(t)
    out = TK.zeros(TK.scalartype(t), TK.space(t, 1) ⊗ E ← TK.domain(t)[1])
    for (f1, f2) in TK.fusiontrees(t)
        if f1.uncoupled[2] == one(G)
            # 物理 (0,0) 块（vac）写进站 (0,0) 槽 1
            out[(f1.uncoupled[1], one(G), TK.dual(only(f2.uncoupled)))][:, 1, :] .= t[f1, f2][:, 1, :]
        else
            # 生成元块写进站形的目标味标签（sector + slot）
            out[(f1.uncoupled[1], tag.sector, TK.dual(only(f2.uncoupled)))][:, tag.slot, :] .= t[f1, f2][:, 1, :]
        end
    end
    out
end

"相邻两站（同味变量对）楔积融合：**收缩中间键** + wedge_vertex 站级楔积顶点；
机器符号由 GM.@grassmann 承担"
function _fuse_adjacent(a, b)
    GM.domain(a)[1] == GM.space(b, 1) ||
        throw(ArgumentError("adjacent-site fusion requires matching middle bond"))
    mv = GM.wedge_vertex(TK.scalartype(a), GM.space(a, 2))
    GM.@grassmann out[l p; r] := a[l sa; m] * b[m sb; r] * mv[p; sa sb]
    out
end

"""
    pairize(spec, tensors) -> GMPS

变量粒度裸张量链（偶数长）→ pair 布局 GMPS：逐对 (θ站, θ̄站) 嵌入后相邻融合
（中间键收缩 + 楔积）。链长减半、物理维 2×4。
"""
function pairize(spec::AbstractSymmetrySpec, tensors::AbstractVector)
    iseven(length(tensors)) || error("pairize requires an even variable count (θ/θ̄ paired): got $(length(tensors))")
    interp = [_convert_site(spec, t, iseven(i)) for (i, t) in enumerate(tensors)]
    data = [_fuse_adjacent(interp[2j-1], interp[2j]) for j in 1:(length(tensors) ÷ 2)]
    GM.GMPS(data)
end

# ---------------------------------------------------------------------------
# 链构造 API（from_core/split_cell/monomial；全部停留在变量粒度层、输出裸
# 张量——pairize 后才进 GrassmannMPS。U1GTEMPO2 同名函数的对称性泛型版）
# ---------------------------------------------------------------------------

function from_core(spec::AbstractSymmetrySpec, data, leftlabels, rightlabels, c)
    t = TK.zeros(ComplexF64, qspace(spec, leftlabels) ⊗ pspace(spec, c) ← qspace(spec, rightlabels))
    z = zero_charge(spec)
    for q in unique(leftlabels), (si, s) in ((1, z), (2, c))
        r = qadd(q, s)
        r in rightlabels || continue
        li = findall(==(q), leftlabels); ri = findall(==(r), rightlabels)
        t[(sector(spec, q), sector(spec, s), TK.dual(sector(spec, r)))] .= data[li, si:si, ri]
    end
    t
end

# Split only fixed-size unit cells, never the full temporal polynomial.
function split_cell(spec::AbstractSymmetrySpec, data::AbstractArray{ComplexF64,3}, leftlabels, rightlabels, charges)
    work = Array(data); labels = leftlabels
    z = zero_charge(spec)
    T0 = typeof(TK.zeros(ComplexF64, qspace(spec, [z]) ⊗ pspace(spec, charges[1]) ← qspace(spec, [z])))
    tensors = T0[]
    for c in charges[1:end-1]
        rest = size(work, 2) ÷ 2; mat = reshape(work, 2length(labels), :)
        rowqs = [qadd(q, s == z ? z : c) for s in (z, c) for q in labels]
        newqs = typeof(z)[]; pieces = []
        for q in sort(unique(rowqs))
            rows = findall(==(q), rowqs)
            Uc, S, V = MAK.svd_compact(mat[rows, :])
            # No numerical singular-value cutoff on factor construction.
            keep = count(>(0), diag(S)); keep == 0 && continue
            push!(pieces, (rows, Uc[:, 1:keep], S[1:keep, 1:keep] * V[1:keep, :]))
            append!(newqs, fill(q, keep))
        end
        isempty(newqs) && throw(ArgumentError("zero cell"))
        core = zeros(ComplexF64, length(labels), 2, length(newqs)); next = zeros(ComplexF64, length(newqs), rest, length(rightlabels)); offset = 0
        for (rows, Uc, R) in pieces
            k = size(Uc, 2); dst = offset+1:offset+k
            reshape(core, 2length(labels), :)[rows, dst] = Uc
            next[dst, :, :] = reshape(R, k, rest, length(rightlabels)); offset += k
        end
        push!(tensors, from_core(spec, core, labels, newqs, c)); work = next; labels = newqs
    end
    push!(tensors, from_core(spec, work, labels, rightlabels, charges[end])); tensors
end

"""
    monomial(spec, charges, positions=Int[]; coeff=1) -> GMPS

单项式链：`positions` 上的变量带各自生成元荷、其余真空；逆序号吸收进 coeff。
空 charges 抛错（零时刻轮廓走 solve 的零时路径）。
"""
function monomial(spec::AbstractSymmetrySpec, charges::AbstractVector, positions::AbstractVector{<:Integer}=Int[]; coeff::Number=1)
    isempty(charges) && throw(ArgumentError("use the zero-time path for empty contours"))
    all(i -> 1 <= i <= length(charges), positions) || throw(BoundsError(charges, positions))
    length(unique(positions)) == length(positions) || throw(ArgumentError("repeated Grassmann variable"))
    all(isfinite, (coeff,)) || throw(ArgumentError("nonfinite coefficient"))
    inv = sum(positions[i] > positions[j] for i in eachindex(positions) for j in i+1:length(positions); init=0)
    z = zero_charge(spec)
    q = z
    data = map(eachindex(charges)) do i
        s = i in positions ? charges[i] : z; r = qadd(q, s)
        t = TK.zeros(ComplexF64, qspace(spec, [q]) ⊗ pspace(spec, charges[i]) ← qspace(spec, [r]))
        t[(sector(spec, q), sector(spec, s), TK.dual(sector(spec, r)))] .= i == 1 ? coeff * (-1)^inv : 1
        q = r; t
    end
    pairize(spec, data)
end

# ---------------------------------------------------------------------------
# 系数提取（站基矢指标 0..3 经 grassmann_basis 定位，不硬编码索引；
# 占据映射：变量序 (2k−1,2k) → 站 k 的 η/η̄ 槽）
# ---------------------------------------------------------------------------

function _occupation_index(p, o1::Int, o2::Int)
    b = GM.grassmann_basis(p)
    tag = o1 == 0 ? (o2 == 0 ? b.vac : b.etabar) : (o2 == 0 ? b.eta : b.pair)
    basis = [(s, j) for s in collect(TK.sectors(p)) for j in 1:TK.dim(p, s)]
    findfirst(isequal((tag.sector, tag.slot)), basis) - 1
end

# 系数抽查直接消费 GM.coefficient_scaled/coefficient（站级 digit 约定：每站一个
# 0..3 基指标，vac0/pair1/eta2/etabar3）；变量级 (o1,o2)→digit 的映射助手留作
# 测量层（Phase 5）私用。
function _var_coefficient_scaled(x::GM.GMPS, occupations)
    L = length(x)
    length(occupations) == 2L || throw(DimensionMismatch("one binary occupation per variable (2 per station)"))
    all(s -> s in (0, 1), occupations) || throw(ArgumentError("binary occupations required"))
    occ = [_occupation_index(GM.physpace(x, k), occupations[2k-1], occupations[2k]) for k in 1:L]
    z = GM.coefficient_scaled(x, occ)
    ScaledComplex(z.logabs, ComplexF64(z.phase))
end

# ---------------------------------------------------------------------------
# 结构查询（旧导出 API 保名）
# ---------------------------------------------------------------------------

bond_dimensions(x::GM.GMPS) = GM.bonddims(x)
charge_dimensions(x::GM.GMPS) = [Dict(charge_label(chain_spec(x), q) => TK.dim(s, q) for q in TK.sectors(s)) for s in GM.bondspaces(x)]
check_charge(x::GM.GMPS) = _check_structure(x)

# 旧缩放字段皮肤：scaling = exp(logscale)（旧 dense = 张量×scaling^L ≡ GM 系数
# = exp(logscale)×张量积）；setscaling! 对应换算
scaling(x::GM.GMPS) = exp(GM.logscale(x))
setscaling!(x::GM.GMPS, s::Number) = GM.setlogscale!(x, log(s))

# ---------------------------------------------------------------------------
# 收缩引擎接线（U1GTEMPO2 §2.5：integrate_scaled→Zvalue、EnvironmentCache→
# environments、中切口分母条件数在本层）
# ---------------------------------------------------------------------------

function integrate_scaled(states::GM.GMPS...)
    1 <= length(states) <= 3 || throw(ArgumentError("integrate supports one to three independent factors"))
    z = GM.Zvalue(collect(states))
    ScaledComplex(z.logabs, ComplexF64(z.phase))
end
integrate(states::GM.GMPS...) = unscale(integrate_scaled(states...))

"""
快照 + GrassmannMPS ExpectationCache 包装 + 中切口分母条件数（旧 middle-cut
cross 语义：|boundary_scalar(L_mid·R_mid)|/norms，站粒度切割取链半长）。
"""
struct EnvironmentCache
    cache::GM.ExpectationCache
    denominator_condition::Float64
end
function EnvironmentCache(states::Tuple)
    1 <= length(states) <= 3 || throw(ArgumentError("one to three factors required"))
    cache = GM.environments(states...)
    mid = cache.L ÷ 2                      # 站粒度半链切割（旧对粒度 mid 同语义）
    lmid = cache.leftenv[mid+1]; rmid = cache.rightenv[mid+1]
    cond = abs(GM._boundary_pair(lmid, rmid)) / max(norm(lmid) * norm(rmid), floatmin(Float64))
    EnvironmentCache(cache, min(cond, cache.denominator_condition))
end
integrate_scaled(c::EnvironmentCache) = integrate_scaled(c.cache.states...)
integrate(c::EnvironmentCache) = unscale(integrate_scaled(c))
